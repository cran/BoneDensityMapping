## To-do
# error message (or fix) for bone model on border of scan?
# does plot cross section need all the bits in example?
# examples in surface normal intersect and voxel point intersect for single bone
# make maxi and mini defaults max and min of vector
# voxel point interesect, add surface mesh

#' import landmark coordinates
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param landmark_path String. File path to landmark data. .json or .fcsv
#' format
#' @return dataframe. Columns are landmark name, x, y, and z coordinates
#' @examples
#' landmark_path <- system.file("extdata", "test_femur.mrk.json",
#'                              package = "BoneDensityMapping")
#' landmarks <- import_lmks(landmark_path)
#' @importFrom rjson fromJSON
#' @importFrom tools file_ext
#' @export
import_lmks <- function(landmark_path) {
  file_type <- file_ext(landmark_path)

  if (file_type == "json") {
    lmks <- fromJSON(file = landmark_path)

    # extract point lists
    lmks_ <- lmks$markups[[1]]$controlPoints

    # extract names and positions
    lmk_names <- rep(NA, length(lmks_))
    coords <- matrix(NA, nrow = length(lmks_), ncol = 3)
    for (i in seq_along(lmks_)) {
      fid <- lmks_[[i]]
      lmk_names[i] <- fid$label
      coords[i, ] <- fid$position
    }
    df <- data.frame(lmk_id = lmk_names, x = coords[, 1], y = coords[, 2], z = coords[, 3])

  } else if (file_type == "fcsv") {
    coords <- read.csv(landmark_path, skip = 3, header = FALSE)[, 2:4]
    df <- data.frame(lmk_id = seq_len(nrow(coords)), x = coords[[1]], y = coords[[2]], z = coords[[3]])

  } else {
    stop("Unsupported file type: must be .json or .fcsv")
  }
  return(df)
}


#' import CT scan
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param scan_filepath String. File path to CT scan data. Should be .nii or .nrrd
#' @return scan object
#' @examples
#' \donttest{
#'   # Download CT scan
#'   url <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/test_CT_hip.nii.gz"
#'   scan_filepath <- tempfile(fileext = ".nii.gz")
#'   download.file(url, scan_filepath, mode = "wb")
#'   import_scan(scan_filepath)
#' }
#' @importFrom oro.nifti readNIfTI nifti
#' @importFrom nat read.nrrd
#' @export
import_scan <- function(scan_filepath) {
  file_type <- tools::file_ext(scan_filepath)

  if (file_type == "nii" | file_type == "gz") {
    nifti_scan <- oro.nifti::readNIfTI(scan_filepath, reorient = FALSE)

  } else if (file_type == "nrrd") {
    nrrd <- read.nrrd(scan_filepath)

    header <- attr(nrrd, "header")
    space_dirs <- header$`space directions`
    origin <- header$`space origin`

    if (is.null(space_dirs) || is.null(origin)) {
      stop("Missing 'space directions' or 'space origin' in NRRD header.")
    }

    # Flip X (L → R) and Y (P → A)
    space_dirs[ , 1:2] <- -space_dirs[ , 1:2]  # flip X and Y axes
    origin[1:2] <- -origin[1:2]                # flip X and Y offsets

    # Convert to NIfTI
    nifti_scan <- oro.nifti::nifti(img = nrrd)

    # Assign affine matrix (srow_x/y/z)
    affine <- matrix(0, nrow = 3, ncol = 4)
    affine[, 1:3] <- space_dirs
    affine[, 4] <- origin

    nifti_scan@srow_x <- affine[1, ]
    nifti_scan@srow_y <- affine[2, ]
    nifti_scan@srow_z <- affine[3, ]

    # Assign qoffsets
    nifti_scan@qoffset_x <- origin[1]
    nifti_scan@qoffset_y <- origin[2]
    nifti_scan@qoffset_z <- origin[3]

    # Tell NIfTI to use affine
    nifti_scan@qform_code <- 1
  } else {
    stop("Unsupported file type: must be .nii or .nrrd")
  }

  return(nifti_scan)
}


#' import surface mesh
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param surface_mesh_filepath String. File path to bone models. .stl or .ply
#' @return mesh object
#' @examples
#' \donttest{
#'   # Download bone model
#'   url <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/test_CT_femur.stl"
#'   bone_filepath <- tempfile(fileext = ".stl")
#'   download.file(url, bone_filepath, mode = "wb")
#'   import_mesh(bone_filepath)
#' }
#' @importFrom Rvcg vcgImport
#' @export
import_mesh <- function(surface_mesh_filepath) {
  # import scan
  surface_mesh <- vcgImport(surface_mesh_filepath)

  return(surface_mesh)
}


#' Check landmarks are close to the mesh
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param surface_mesh mesh object
#' @param landmarks Dataframe. Columns are landmark name, x, y, and z coords
#' @param threshold Numeric. Distance landmark can be from surface without
#' warning being thrown
#' @return String. Returns a message warning that landmarks are not on bone
#' surface
#' @examples
#' \donttest{
#'   # Download bone model
#'   url <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/test_CT_femur.stl"
#'   bone_filepath <- tempfile(fileext = ".stl")
#'   download.file(url, bone_filepath, mode = "wb")
#'   surface_mesh <- import_mesh(bone_filepath)
#'   landmark_path <- system.file("extdata", "test_femur.fcsv",
#'                              package = "BoneDensityMapping")
#'   landmarks <- import_lmks(landmark_path)
#'   landmark_check(surface_mesh, landmarks, threshold = 1.0)
#' }
#' @importFrom rdist cdist
#' @importFrom utils read.csv
#' @export
landmark_check <- function(surface_mesh, landmarks, threshold = 1.0) {
  vertices <- t(surface_mesh$vb)[, c(1:3)]
  coords <- landmarks[, c("x", "y", "z")]

  dists <- c()
  for (i in 1:nrow(coords)) {
    pt <- matrix(as.numeric(unlist(coords[i, ])), nrow = 1)
    x <- cdist(vertices, pt)
    dists <- c(dists, min(x))
  }

  if (any(dists > threshold)) {
    bad_ids <- landmarks$lmk_id[dists > threshold]
    message("Landmarks not on bone surface: ", paste(bad_ids, collapse = ", "))
  }
  else message("All landmarks are on bone surface.")
}


#' Check if surface model is fully contained within scan volume
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param surface_mesh mesh object (class \code{mesh3d}) or numeric
#'    matrix/dataframe of vertex coordinates (cols: X, Y, Z)
#' @param nifti NIfTI image object representing CT scan.
#' @param return_limits Logical. If TRUE returns a summary of the bounding boxes
#'  of the scan and mesh
#' @return If any vertices lie outside the scan volume, it
#'   raises an error.
#' @examples
#' \donttest{
#'   # Download CT scan
#'   url <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/test_CT_hip.nii.gz"
#'   scan_filepath <- tempfile(fileext = ".nii.gz")
#'   download.file(url, scan_filepath, mode = "wb")
#'   nifti <- import_scan(scan_filepath)
#'   url2 <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/test_CT_femur.stl"
#'   bone_filepath <- tempfile(fileext = ".stl")
#'   download.file(url2, bone_filepath, mode = "wb")
#'   surface_mesh <- import_mesh(bone_filepath)
#'   bone_scan_check(surface_mesh, nifti, return_limits = TRUE)
#' }
#' @export
bone_scan_check <- function(surface_mesh, nifti, return_limits = FALSE) {
  # check input mesh
  if (inherits(surface_mesh, "mesh3d")) {
    vertices <- t(surface_mesh$vb)[, 1:3]
  } else if (is.matrix(surface_mesh) || is.data.frame(surface_mesh)) {
    vertices <- as.matrix(surface_mesh)
  } else {
    stop("surface_mesh must be a mesh3d object or a matrix of vertex coordinates.")
  }

  # format image data, with voxel coordinates
  img_data <- img_data(nifti)
  dims <- dim(img_data)
  x_seq <- seq((nifti)@qoffset_x * -1, by = (nifti)@srow_x[1] * -1, length.out = dims[1])
  y_seq <- seq((nifti)@qoffset_y * -1, by = (nifti)@srow_y[2] * -1, length.out = dims[2])
  z_seq <- seq((nifti)@qoffset_z, by = (nifti)@srow_z[3], length.out = dims[3])

  # voxel sizes
  voxel_x <- abs(nifti@srow_x[1])
  voxel_y <- abs(nifti@srow_y[2])
  voxel_z <- abs(nifti@srow_z[3])

  # mesh bounds
  mesh_x_min <- min(vertices[, 1])
  mesh_x_max <- max(vertices[, 1])
  mesh_y_min <- min(vertices[, 2])
  mesh_y_max <- max(vertices[, 2])
  mesh_z_min <- min(vertices[, 3])
  mesh_z_max <- max(vertices[, 3])

  # scan bounds
  vol_x_min <- min(x_seq)
  vol_x_max <- max(x_seq)
  vol_y_min <- min(y_seq)
  vol_y_max <- max(y_seq)
  vol_z_min <- min(z_seq)
  vol_z_max <- max(z_seq)

  # check for outside bounds
  x1_diff <- vol_x_min - mesh_x_min
  x2_diff <- mesh_x_max - vol_x_max
  y1_diff <- vol_y_min - mesh_y_min
  y2_diff <- mesh_y_max - vol_y_max
  z1_diff <- vol_z_min - mesh_z_min
  z2_diff <- mesh_z_max - vol_z_max

  diffs <- c(x1_diff, x2_diff, y1_diff, y2_diff, z1_diff, z2_diff)
  voxel_sizes <- c(voxel_x, voxel_x, voxel_y, voxel_y, voxel_z, voxel_z)

  outside_flags <- diffs > 0
  exceeds_voxel <- diffs > voxel_sizes

  if (any(outside_flags)) {
    if (any(exceeds_voxel)) {
      stop("Mesh not within scan volume.")
    } else {
      message("Mesh is within 1 voxel outside scan volume. Consider cropping mesh accordingly.")
    }
  }

  df <- data.frame(
    Axis = c("X", "Y", "Z"),
    Mesh_Min = c(mesh_x_min, mesh_y_min, mesh_z_min),
    Mesh_Max = c(mesh_x_max, mesh_y_max, mesh_z_max),
    Scan_Min = c(vol_x_min, vol_y_min, vol_z_min),
    Scan_Max = c(vol_x_max, vol_y_max, vol_z_max)
  )

  if (return_limits == TRUE) {
    return(df)
  }
}



#' Fills bone with orthogonally spaced points for internal analysis
#' @param surface_mesh Mesh object
#' @param spacing Numeric
#' @return Matrix with internal point coordinates
#' @examples
#' \donttest{
#'   url <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/test_CT_femur.stl"
#'   bone_filepath <- tempfile(fileext = ".stl")
#'   download.file(url, bone_filepath, mode = "wb")
#'   surface_mesh <- import_mesh(bone_filepath)
#'   internal_fill <- fill_bone_points(surface_mesh, 2)
#' }
#' @importFrom ptinpoly pip3d
#' @importFrom stats runif
#' @export
fill_bone_points <- function(surface_mesh, spacing) {
  vertices <- t(surface_mesh$vb)[, 1:3]
  faces <- t(surface_mesh$it)

  # bone extents
  x_min <- min(vertices[, 1])
  x_max <- max(vertices[, 1])
  y_min <- min(vertices[, 2])
  y_max <- max(vertices[, 2])
  z_min <- min(vertices[, 3])
  z_max <- max(vertices[, 3])

  # make point df
  x <- seq(from = x_min, to = x_max, by = spacing)
  y <- seq(from = y_min, to = y_max, by = spacing)
  z <- seq(from = z_min, to = z_max, by = spacing)
  pt_mat <- as.matrix(expand.grid(x = x, y = y, z = z))

  # find which points are within surface
  in_bone <- which(pip3d(vertices, faces, pt_mat) == 1)
  in_coords <- pt_mat[in_bone, ]

  # add a little noise
  dims <- dim(in_coords)
  noise <- runif(dims[1] * dims[2], -0.001, 0.001)
  noise_mat <- matrix(noise, nrow = dims[1], ncol = dims[2])
  in_coords <- in_coords + noise_mat

  return(in_coords)
}


#' Redefine surface points. Adds additional surface points (“sliders”) that
#' are spatially distributed across the mesh surface. Adapted from geomorph
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param surface_mesh Mesh object
#' @param landmarks Data frame with landmark coordinates (columns: ID, x, y, z)
#' @param no_surface_sliders Numeric. No. of surface points to generate
#' @return Data frame. 3D coordinates for the combined set of original
#' landmarks and the new surface points
#' @examples
#' \donttest{
#'   url <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/test_CT_femur.stl"
#'   bone_filepath <- tempfile(fileext = ".stl")
#'   download.file(url, bone_filepath, mode = "wb")
#'   surface_mesh <- import_mesh(bone_filepath)
#'   landmark_path <- system.file("extdata", "test_femur.mrk.json",
#'                                package = "BoneDensityMapping")
#'   landmarks <- import_lmks(landmark_path)
#'   mapped_coords <- surface_points_template(surface_mesh, landmarks, 1000)
#' }
#' @importFrom stats kmeans
#' @export
surface_points_template <- function(surface_mesh, landmarks,
                                    no_surface_sliders) {
  # extract vertex coordinates from mesh
  if (is.list(surface_mesh)) {
    vertices <- t(surface_mesh$vb)[, c(1:3)]
  } else {
    vertices <- surface_mesh
  }
  colnames(vertices) <- c("xpts", "ypts", "zpts")

  # Isolate just the x, y, z coordinates from the landmark df
  landmark_coords <- landmarks[, 2:4]
  landmark_coords <- apply(landmark_coords, 2, as.numeric)  # ensure numeric

  # Identify which mesh vertices are closest to each landmark, add them to lmk_add, remove these mesh vertices
  lmk.add <- NULL
  for(i in 1:nrow(landmark_coords)){
    lmk.add <- rbind(lmk.add, which.min(sqrt((landmark_coords[i, 1] - vertices[, 1]) ^ 2 +
                                               (landmark_coords[i, 2] - vertices[, 2]) ^ 2 +
                                               (landmark_coords[i, 3] - vertices[, 3]) ^ 2))[1])
  }
  vertices <- vertices[-lmk.add, ]

  # Use k-means clustering to find 'no_surface_sliders' new surface points
  colnames(landmark_coords) <- c("xpts", "ypts", "zpts")
  new_surface_points <- rbind(landmark_coords,
                              kmeans(x = vertices, centers = no_surface_sliders,
                                     iter.max = 25)$centers)

  return(as.data.frame(new_surface_points))
}


#' New mapped surface points from template
#' @author Scott Telfer \email{scott.telfer@gmail.com} Adapted from geomorph
#' @param surface_mesh List. Mesh data imported via ply_import function
#' @param landmarks Data frame. Contains 3D coords of landmarks
#' @param template Data frame. 3D coords of remapped surface points
#' @param mirror Logical or character. Set to "x", "y", or "z" to mirror the
#'  mesh and landmarks across that axis before remapping.
#' @param plot_check Logical. If TRUE, generates a 3D plot showing the mirrored
#'  mesh, mirrored landmarks, remapped surface points, and original template
#'  points to visually verify correct orientation and laterality.
#' @return Data frame. 3D coords of remapped surface points
#' @examples
#' \donttest{
#'   url <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/SCAP001.stl"
#'   bone_filepath <- tempfile(fileext = ".stl")
#'   download.file(url, bone_filepath, mode = "wb")
#'   scap_001_mesh <- import_mesh(bone_filepath)
#'   landmark_path <- system.file("extdata", "SCAP001_landmarks.fcsv",
#'                                package = "BoneDensityMapping")
#'   scap_001_lmk <- import_lmks(landmark_path)
#'   template_coords <- surface_points_template(scap_001_mesh, scap_001_lmk,
#'                                              1000)
#'   url2 <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/SCAP002.stl"
#'   bone_filepath <- tempfile(fileext = ".stl")
#'   download.file(url2, bone_filepath, mode = "wb")
#'   scap_002_mesh <- import_mesh(bone_filepath)
#'   landmark_path <- system.file("extdata", "SCAP002_landmarks.fcsv",
#'                                package = "BoneDensityMapping")
#'   scap_002_lmk <- import_lmks(landmark_path)
#'   scap_002_remapped <- surface_points_new(scap_002_mesh, scap_002_lmk,
#'                                           template_coords, mirror = "x",
#'                                           plot_check = FALSE)
#' }
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom rgl open3d points3d title3d
#' @importFrom stats dist
#' @export
surface_points_new <- function(surface_mesh, landmarks, template,
                               mirror = FALSE, plot_check = FALSE) {
  # helper functions
  rotate.mat <- function(M, Y){
    k <- ncol(M)
    M <- cs.scale(M); Y <- cs.scale(Y)
    MY <- crossprod(M, Y)
    sv <- La.svd(MY, k, k)
    u <- sv$u; u[, k] <- u[, k] * determinant(MY)$sign
    v <- t(sv$vt)
    tcrossprod(v, u)
  }

  csize <- function(x) sqrt(sum(center(as.matrix(x))^2))

  center <- function(x){
    if(is.vector(x)) x - mean(x) else {
      x <- as.matrix(x)
      dims <- dim(x)
      fast.center(x, dims[1], dims[2])
    }
  }

  fast.center <- function(x, n, p){
    m <- colMeans(x)
    x - rep.int(m, rep_len(n, p))
  }

  cs.scale <- function(x) x/csize(x)

  fast.solve <- function(x) {
    x <- as.matrix(x)
    if(det(x) > 1e-8) {
      res <- try(chol2inv(chol(x)), silent = TRUE)
      if (inherits(res, "try-error")) res <- fast.ginv(x)
    } else res <- fast.ginv(x)
    return(res)
  }

  tps2d3d <- function(M, matr, matt, PB = TRUE){		#DCA: altered from J. Claude 2008
    p <- dim(matr)[1]; k <- dim(matr)[2]; q <- dim(M)[1]
    Pdist <- as.matrix(stats::dist(matr))
    ifelse(k == 2, P <- Pdist^2*log(Pdist^2), P <- Pdist)
    P[which(is.na(P))] <- 0
    Q <- cbind(1, matr)
    L <- rbind(cbind(P, Q), cbind(t(Q), matrix(0, k + 1, k + 1)))
    m2 <- rbind(matt, matrix(0, k + 1, k))
    coefx <- fast.solve(L)%*%m2[, 1]
    coefy <- fast.solve(L)%*%m2[, 2]
    if(k == 3){coefz <- fast.solve(L)%*%m2[, 3]}
    fx <- function(matr, M, coef, step){
      Xn <- numeric(q)
      for (i in 1:q){
        Z <- apply((matr-matrix(M[i,], p, k, byrow = TRUE))^2, 1, sum)
        ifelse(k == 2, Z1<-Z*log(Z), Z1<-sqrt(Z)); Z1[which(is.na(Z1))] <- 0
        ifelse(k == 2, Xn[i] <- coef[p+1] + coef[p+2]*M[i,1] + coef[p+3]*M[i,2] + sum(coef[1:p]*Z1),
               Xn[i] <- coef[p+1] + coef[p+2]*M[i,1] + coef[p+3]*M[i,2] + coef[p+4]*M[i,3] + sum(coef[1:p]*Z1))
        if(PB == TRUE){setTxtProgressBar(pb, step + i)}
      }
      return(Xn)
    }
    matg <- matrix(NA, q, k)
    if(PB==TRUE){pb <- txtProgressBar(min = 0, max = q*k, style = 3) }
    matg[,1] <- fx(matr, M, coefx, step = 1)
    matg[,2] <- fx(matr, M, coefy, step=q)
    if(k==3){matg[,3] <- fx(matr, M, coefz, step=q*2)
    }
    if(PB==TRUE) close(pb)
    return(matg)
  }

  fast.ginv <- function(X, tol = sqrt(.Machine$double.eps)){
    X <- as.matrix(X)
    k <- ncol(X)
    Xsvd <- La.svd(X, k, k)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
    rtu <-((1 / Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE]))
    v <-t(Xsvd$vt)[, Positive, drop = FALSE]
    v%*%rtu
  }

  # copy inputs to avoid modifying them
  mesh_input <- surface_mesh
  lmk_input <- landmarks

  # apply mirroring
  if (mirror %in% c("x", "y", "z")) {
    axis_idx <- switch(mirror,
                       x = 1,
                       y = 2,
                       z = 3)

    # Mirror mesh3d or matrix mesh
    if (is.list(mesh_input) && !is.null(mesh_input$vb)) {
      mesh_input$vb[axis_idx, ] <- -mesh_input$vb[axis_idx, ]
    } else {
      mesh_input[, axis_idx] <- -mesh_input[, axis_idx]
    }

    # Mirror landmarks dataframe columns x/y/z accordingly
    axis_name <- mirror
    lmk_input[[axis_name]] <- -lmk_input[[axis_name]]
  }

  # format mesh
  if (is.list(mesh_input)) {
    bone <- t(mesh_input$vb)[, c(1:3)]
  } else {
    bone <- mesh_input
  }

  # closest vertex to landmark
  lmk.add <- NULL
  landmark_xyz <- as.matrix(lmk_input[, c("x", "y", "z")])
  for(i in 1:nrow(lmk_input)){
    lmk.add <- rbind(lmk.add,
                     which.min(sqrt((landmark_xyz[i, 1] - bone[, 1]) ^ 2 +
                                      (landmark_xyz[i, 2] - bone[, 2]) ^ 2 +
                                      (landmark_xyz[i, 3] - bone[, 3]) ^ 2))[1])
  }
  nlandmarks <- nrow(lmk_input)

  # center bone
  bone_centered <- center(bone)
  bone_trans <- colMeans(bone)

  # center template
  template <- center(template) * (csize(bone_centered[lmk.add, ]) / csize(template[(1:nlandmarks), ]))
  template <- template %*% rotate.mat(bone_centered[lmk.add, ], template[(1:nlandmarks), ])

  # sliding points
  template.tps <- tps2d3d(template[-(1:nlandmarks), ], template[(1:nlandmarks), ], bone_centered[lmk.add, ])
  spec.surfs <- bone_centered[-lmk.add, ]
  nei <- numeric(dim(template.tps)[1])
  sliders <- matrix(NA, nrow = dim(template.tps)[1], ncol = 3)
  for (i in 1:dim(template.tps)[1]) {
    nei[i] <- which.min(sqrt((template.tps[i, 1] - spec.surfs[, 1]) ^ 2 +
                               (template.tps[i, 2] - spec.surfs[ ,2]) ^ 2 +
                               (template.tps[i, 3] - spec.surfs[, 3]) ^ 2))[1]
    sliders[i,] <- spec.surfs[nei[i], ]
    spec.surfs <- spec.surfs[-nei[i], ]
  }

  # make output matrix
  selected.out <- rbind(bone_centered[lmk.add, ], sliders)

  # translate back
  selected.out[, 1] <- selected.out[, 1] - (bone_trans[1] * - 1)
  selected.out[, 2] <- selected.out[, 2] - (bone_trans[2] * - 1)
  selected.out[, 3] <- selected.out[, 3] - (bone_trans[3] * - 1)

  if (!identical(mirror, FALSE) && plot_check) {
    open3d()
    shade3d(mesh_input, color = "gray", alpha = 0.2)
    points3d(as.matrix(lmk_input[, c("x", "y", "z")]), col = "red", size = 2)
    points3d(selected.out, col = "purple", size = 0.5)
    points3d(template, col = "green", size = 0.5)
    title3d(main = "Check: Laterality Alignment Before Un-Mirroring // green = template", col = "black", line = 2)
  }

  # unmirror if needed
  if (mirror %in% c("x", "y", "z")) {
    axis_idx <- switch(mirror, x = 1, y = 2, z = 3)
    selected.out[, axis_idx] <- -selected.out[, axis_idx]
  }

  if (plot_check) {
    open3d()
    shade3d(surface_mesh, color = "gray", alpha = 0.2)
    points3d(as.matrix(landmarks[, c("x", "y", "z")]), col = "red", size = 2)
    points3d(selected.out, col = "purple", size = 0.5)
    points3d(template, col = "green", size = 0.5)
    title3d(main = "Check: Completed new surface points", col = "black", line = 2)

  }

  return(selected.out)
}


#' Find material properties of bone at surface point using surface normal
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param surface_mesh Mesh object
#' @param mapped_coords Data frame. 3D coords of remapped surface points. If
#' NULL, surface_mesh vertices will be used
#' @param normal_dist Numeric. Distance surface normal should penetrate surface
#' @param nifti Nifti CT scan image
#' @param ct_eqn String. Equation to use for density calibration. Currently
#' "linear" supported.
#' @param ct_params Numeric vector. Calibration parameters for density
#' calculation. For linear, first value is beta coefficient (y intercept),
#' second value is sigma coefficient (gradient)
#' @param rev_x Logical. Reverses x voxel coordinates
#' @param rev_y Logical. Reverses y voxel coordinates
#' @param rev_z Logical. Reverses z voxel coordinates
#' @param check_in_vol Logical. Include check that model is in scans volume
#' and print dimensions
#' @return Vector. Vector with value for each point on surface
#' @examples
#' \donttest{
#'   # Download CT scan
#'   url <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/test_CT_hip.nii.gz"
#'   scan_filepath <- tempfile(fileext = ".nii.gz")
#'   download.file(url, scan_filepath, mode = "wb")
#'   nifti <- import_scan(scan_filepath)
#'   url2 <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/test_CT_femur.stl"
#'   bone_filepath <- tempfile(fileext = ".stl")
#'   download.file(url2, bone_filepath, mode = "wb")
#'   surface_mesh <- import_mesh(bone_filepath)
#'   landmark_path <- system.file("extdata", "test_femur.mrk.json",
#'                                package = "BoneDensityMapping")
#'   landmarks <- import_lmks(landmark_path)
#'   mapped_coords <- surface_points_template(surface_mesh, landmarks,
#'                                            no_surface_sliders = 1000)
#'   mat_peak <- surface_normal_intersect(surface_mesh, normal_dist = 3.0,
#'                                        nifti = nifti, ct_eqn = "linear",
#'                                        ct_params = c(68.4, 1.106))
#' }
#' @importFrom oro.nifti img_data
#' @importFrom RNifti niftiHeader
#' @export
surface_normal_intersect <- function(surface_mesh, mapped_coords = NULL,
                                     normal_dist = 3.0, nifti,
                                     ct_eqn = NULL, ct_params = NULL,
                                     rev_x = FALSE, rev_y = FALSE,
                                     rev_z = FALSE, check_in_vol = FALSE) {
  # extract details from mesh
  surface_coords <- t(surface_mesh$vb)[, c(1:3)]
  surface_normals <- t(surface_mesh$normals)[, c(1:3)]

  # Convert mapped_coords to numeric matrix for calculations
  if (hasArg(mapped_coords)) {
    vertices <- data.matrix(mapped_coords)
    dims <- dim(vertices)
    vertices <- as.numeric(vertices)
    dim(vertices) <- dims
    mat_peak <- rep(NA, times = nrow(vertices))
  } else {
    mat_peak <- rep(NA, times = nrow(surface_coords))
    vertices <- surface_coords
  }

  # format image data, with voxel coordinates
  img_data <- img_data(nifti)
  dims <- dim(img_data)
  x_by <- (nifti)@srow_x[1]
  y_by <- (nifti)@srow_y[2]
  z_by <- (nifti)@srow_z[3]
  if (rev_x == TRUE) {
    x_seq <- rev(seq(nifti@qoffset_x * -1, by = x_by * -1, length.out = dims[1]))
  } else {
    x_seq <- seq((nifti)@qoffset_x * -1, by = x_by * -1, length.out = dims[1])
  }
  if (rev_y == TRUE) {
    y_seq <- rev(seq((nifti)@qoffset_y * -1, by = y_by * -1, length.out = dims[2]))
  } else {
    y_seq <- seq((nifti)@qoffset_y * -1, by = y_by * -1, length.out = dims[2])
  }
  if (rev_z == TRUE) {
    z_seq <- rev(seq((nifti)@qoffset_z, by = z_by, length.out = dims[3]))
  } else {
    z_seq <- seq((nifti)@qoffset_z, by = z_by, length.out = dims[3])
  }

  # check vertices are in scan volume
  if (check_in_vol == TRUE) {
    bone_scan_check(surface_coords, nifti)
  }

  # Find voxels intercepted by line
  ## Loop through each surface point to sample CT data along the surface normal
  for (i in 1:nrow(vertices)) {
    # find start and end points of normal line
    if (hasArg(mapped_coords)) {
      yy <- t(as.matrix(vertices[i, ], 1, 3))
      y <- cdist(yy, surface_coords)
      matched_point <- which.min(y)
      start_point <- surface_coords[matched_point, ]
      end_point <- surface_coords[matched_point, ] + (surface_normals[matched_point, ] * -1 * normal_dist)
    } else {
      start_point <- surface_coords[i, ]
      end_point <- surface_coords[i, ] + (surface_normals[i, ] * -1 * normal_dist)
    }

    # Interpolate 10 points along this line inside the bone surface
    px <- seq(from = start_point[1], to = end_point[1], length.out = 10)
    py <- seq(from = start_point[2], to = end_point[2], length.out = 10)
    pz <- seq(from = start_point[3], to = end_point[3], length.out = 10)

    # For each point along the line, find corresponding voxel and record CT intensity
    max_line <- rep(NA, times = 10)
    for (j in 1:10) {
      voxel <- c(which.min(abs(px[j] - x_seq)),
                 which.min(abs(py[j] - y_seq)),
                 which.min(abs(pz[j] - z_seq)))

      max_line[j] <- img_data[voxel[1], voxel[2], voxel[3]]
    }

    # add HU column
    mat_peak[i] <- max(max_line)
  }

  # Convert CT numbers to density using calibration values
  if (hasArg(ct_eqn)) {
    mat_peak <- ct_calibration(mat_peak, ct_eqn, ct_params)
  }

  # Return the vector of density values for each mapped coordinate
  return(mat_peak)
}


#' Finds material properties of bone at any point
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param vertex_coords Matrix
#' @param nifti nifti object
#' @param ct_eqn String. Equation to use for density calibration. Currently
#' only "linear" is supported.
#' @param ct_params Numeric vector. Calibration parameters for density
#' calculation. For linear, first value is beta coefficient (y intercept),
#' second value is sigma coefficient (gradient)
#' @param check_in_vol Logical. Include check that model is in scans volume
#' and print dimensions
#' @param rev_x Logical. Reverses x voxel coordinates
#' @param rev_y Logical. Reverses y voxel coordinates
#' @param rev_z Logical. Reverses z voxel coordinates
#' @examples
#' \donttest{
#'   # Download CT scan
#'   url <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/test_CT_hip.nii.gz"
#'   scan_filepath <- tempfile(fileext = ".nii.gz")
#'   download.file(url, scan_filepath, mode = "wb")
#'   nifti <- import_scan(scan_filepath)
#'   url2 <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/test_CT_femur.stl"
#'   bone_filepath <- tempfile(fileext = ".stl")
#'   download.file(url2, bone_filepath, mode = "wb")
#'   surface_mesh <- import_mesh(bone_filepath)
#'
#'   # get density of surface bone directly
#'   mat_peak <- voxel_point_intersect(surface_mesh, nifti,
#'                                     ct_eqn = "linear",
#'                                     ct_params = c(68.4, 1.106),
#'                                     check_in_vol = FALSE)
#'
#'   # remap and get density (for group level comparisons)
#'   landmark_path <- system.file("extdata", "test_femur.mrk.json",
#'                                package = "BoneDensityMapping")
#'   landmarks <- import_lmks(landmark_path)
#'   mapped_coords <- surface_points_template(surface_mesh, landmarks,
#'                                            no_surface_sliders = 1000)
#'   mat_peak <- voxel_point_intersect(mapped_coords, nifti,
#'                                     ct_eqn = "linear",
#'                                     ct_params = c(68.4, 1.106),
#'                                     check_in_vol = FALSE)
#' }
#' @return Vector. Vector with value for each point on surface
#' @export
voxel_point_intersect <- function(vertex_coords, nifti, ct_eqn = NULL,
                                  ct_params = NULL, rev_x = FALSE,
                                  rev_y = FALSE, rev_z = FALSE,
                                  check_in_vol = FALSE) {
  # check input mesh
  if (inherits(vertex_coords, "mesh3d")) {
    vert_coords <- t(vertex_coords$vb)[, 1:3]
  } else if (is.matrix(vertex_coords) || is.data.frame(vertex_coords)) {
    vert_coords <- as.matrix(vertex_coords)
  } else {
    stop("surface_mesh must be a mesh3d object or a matrix of vertex coordinates.")
  }

  #vertex_coords <- data.matrix(vertex_coords)
  #dims <- dim(vertex_coords)
  #vertex_coords <- as.numeric(vertex_coords)
  #dim(vertex_coords) <- dims

  # format image data, with voxel coordinates
  img_data <- img_data(nifti)
  dims <- dim(img_data)
  x_by <- (nifti)@srow_x[1]
  y_by <- (nifti)@srow_y[2]
  z_by <- (nifti)@srow_z[3]
  if (rev_x == TRUE) {
    x_seq <- rev(seq(nifti@qoffset_x * -1, by = x_by * -1, length.out = dims[1]))
  } else {
    x_seq <- seq((nifti)@qoffset_x * -1, by = x_by * -1, length.out = dims[1])
  }
  if (rev_y == TRUE) {
    y_seq <- rev(seq((nifti)@qoffset_y * -1, by = y_by * -1, length.out = dims[2]))
  } else {
    y_seq <- seq((nifti)@qoffset_y * -1, by = y_by * -1, length.out = dims[2])
  }
  if (rev_z == TRUE) {
    z_seq <- rev(seq((nifti)@qoffset_z, by = z_by, length.out = dims[3]))
  } else {
    z_seq <- seq((nifti)@qoffset_z, by = z_by, length.out = dims[3])
  }

  ## check vertices are in scan volume
  if (check_in_vol == TRUE) {
    bone_scan_check(vert_coords, nifti)
  }

  ## Find voxels intercepted by line
  mat_peak <- rep(NA, times = nrow(vert_coords))
  for (i in 1:nrow(vert_coords)) {
    # point coordinates
    point <- vert_coords[i, ]

    # find voxel
    voxel <- c(which.min(abs(point[1] - x_seq)),
               which.min(abs(point[2] - y_seq)),
               which.min(abs(point[3] - z_seq)))

    # add HU column
    mat_peak[i] <- img_data[voxel[1], voxel[2], voxel[3]]
  }

  # Convert CT numbers to density using calibration values
  if (hasArg(ct_eqn)) {
    mat_peak <- ct_calibration(mat_peak, "linear", ct_params)
  }

  # return
  return(mat_peak)
}


#' maps numeric values to a color
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param x Vector.
#' @param maxi Numeric. Maximum value. Defaults to maximum value in vector.
#' Can be useful to set this manually if there are outliers in the scan data
#' @param mini Numeric. Minimum value. Defaults to minimum value in vector.
#' Can be useful to set this manually if there are outliers in the scan data
#' @param color_sel Vector. Colors to use for map. Defaults to a scale of
#' "grey", "blue", "green", "yellow", "orange", "red", "pink".
#' @examples
#' \donttest{
#'   # Download CT scan
#'   url <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/test_CT_hip.nii.gz"
#'   scan_filepath <- tempfile(fileext = ".nii.gz")
#'   download.file(url, scan_filepath, mode = "wb")
#'   nifti <- import_scan(scan_filepath)
#'   url2 <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/test_CT_femur.stl"
#'   bone_filepath <- tempfile(fileext = ".stl")
#'   download.file(url2, bone_filepath, mode = "wb")
#'   surface_mesh <- import_mesh(bone_filepath)
#'   mat_peak <- voxel_point_intersect(surface_mesh, nifti, ct_eqn = "linear",
#'                                     ct_params = c(68.4, 1.106),)
#'   colors <- color_mapping(mat_peak)
#' }
#' @return Vector of hex color values same length as x
#' @importFrom grDevices colorRamp rgb
#' @export
color_mapping <- function(x, maxi, mini, color_sel) {
  # Scale input vector x between 0 and 1
  if (missing(maxi) & missing(mini)) {
    x01 <- (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  } else if (!missing(maxi) & !missing(mini)) {
    x01 <- (x - mini) / (maxi - mini)
  } else {
    stop("Either provide both 'maxi' and 'mini' or neither.")
  }

  # Replace NA values with 0
  x01[is.na(x01)] <- 0

  # Clamp values between 0 and 1 to avoid colorRamp errors
  x01[x01 < 0] <- 0
  x01[x01 > 1] <- 1

  # Generate color map
  if (missing(color_sel)) {
    colormap <- colorRamp(c("grey", "blue", "green", "yellow", "orange", "red", "pink"),
                          interpolate = "spline")
  } else {
    colormap <- colorRamp(color_sel)
  }

  # Apply color map to scaled values
  ply_col <- colormap(x01)
  ply_col <- apply(ply_col, 1, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))

  return(ply_col)
}


#' Takes a density vector mapped to standardized coordinates and maps it to a
#' surface mesh for visualization.
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param surface_mesh Mesh object
#' @param template_pts Matrix
#' @param density_vector Vector
#' @param maxi Numeric
#' @param mini Numeric
#' @param export_path Character
#' @param color_sel String
#' @return mesh3d object with added color dimension
#' @examples
#' \donttest{
#'   # Download CT scan
#'   url <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.1/test_CT_hip.nii.gz"
#'   scan_filepath <- tempfile(fileext = ".nii.gz")
#'   download.file(url, scan_filepath, mode = "wb")
#'   nifti <- import_scan(scan_filepath)
#'   url2 <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/test_CT_femur.stl"
#'   bone_filepath <- tempfile(fileext = ".stl")
#'   download.file(url2, bone_filepath, mode = "wb")
#'   surface_mesh <- import_mesh(bone_filepath)
#'   landmark_path <- system.file("extdata", "test_femur.mrk.json",
#'                                package = "BoneDensityMapping")
#'   landmarks <- import_lmks(landmark_path)
#'   mapped_coords <- surface_points_template(surface_mesh, landmarks,
#'                                            no_surface_sliders = 1000)
#'   dens <- voxel_point_intersect(mapped_coords, nifti,
#'                                 ct_eqn = "linear",
#'                                 ct_params = c(68.4, 1.106))
#'   colored_mesh <- color_mesh(surface_mesh, mapped_coords, dens)
#' }
#' @importFrom Rvcg vcgPlyWrite
#' @importFrom FNN get.knnx
#' @export
color_mesh <- function(surface_mesh, template_pts, density_vector,
                       maxi = NULL, mini = NULL, export_path, color_sel) {

  mesh_template_match <- function(surface_mesh, template_points) {
    # Get vertex coordinates from the mesh
    vertex_coords <- t(surface_mesh$vb)[, c(1:3)]
    template_points <- as.matrix(template_points)

    # nearest neighbor match
    nn <- get.knnx(data = template_points, query = vertex_coords, k = 1)
    return(as.vector(nn$nn.index))
  }

  mesh_match <- mesh_template_match(surface_mesh, template_pts)

  # Pass maxi and mini only if both are not NULL, else call without them
  if (!is.null(maxi) && !is.null(mini)) {
    color_map <- color_mapping(density_vector, maxi, mini, color_sel = color_sel)
  } else {
    color_map <- color_mapping(density_vector, color_sel = color_sel)
  }

  # color to mesh
  surface_mesh$material$color <- color_map[mesh_match]

  # return
  if (missing(export_path) == FALSE) {
    vcgPlyWrite(surface_mesh, export_path)
  } else {
    return(surface_mesh)
  }
}


#' plot mesh
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param surface_mesh Mesh object
#' @param density_color Vector. Colors mapped from density values.
#' @param title String. Plot title.
#' @param userMat Optional matrix. Controls graph orientation.
#' @param legend Logical. Optional color bar.
#' @param legend_color_sel Optional character with color gradient
#' @param legend_maxi Numeric. Maximum bone density.
#' @param legend_mini Numeric. Minimum bone density.
#' @return plot of mesh with color
#' @examples
#' \donttest{
#'   # Download CT scan
#'   url <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.1/test_CT_hip.nii.gz"
#'   scan_filepath <- tempfile(fileext = ".nii.gz")
#'   download.file(url, scan_filepath, mode = "wb")
#'   nifti <- import_scan(scan_filepath)
#'   url2 <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/test_CT_femur.stl"
#'   bone_filepath <- tempfile(fileext = ".stl")
#'   download.file(url2, bone_filepath, mode = "wb")
#'   surface_mesh <- import_mesh(bone_filepath)
#'   vertices <- t(surface_mesh$vb)[, c(1:3)]
#'   landmark_path <- system.file("extdata", "test_femur.mrk.json",
#'                                package = "BoneDensityMapping")
#'   landmarks <- import_lmks(landmark_path)
#'   mapped_coords <- surface_points_template(surface_mesh, landmarks,
#'                                            no_surface_sliders = 5000)
#'   mat_peak <- surface_normal_intersect(surface_mesh, mapped_coords,
#'                                        normal_dist = 3.0, nifti, rev_y=FALSE)
#'   color_mesh <- color_mesh(surface_mesh, mapped_coords, mat_peak, maxi = 2000,
#'                            mini = 0)
#'   plot <- plot_mesh(color_mesh)
#' }
#' @importFrom rgl shade3d view3d bgplot3d
#' @importFrom graphics plot.new mtext
#' @importFrom methods hasArg
#' @importFrom grDevices colorRampPalette
#' @export
plot_mesh <- function(surface_mesh, density_color = NULL, title = NULL,
                      legend = TRUE, legend_color_sel = NULL,
                      legend_maxi = 2000, legend_mini = 0,
                      userMat = NULL) {

  vertices <- as.matrix(t(surface_mesh$vb)[,-4])
  surface_mesh$vb <- rbind(t(vertices), 1)

  open3d()
  par3d(windowRect = c(20, 30, 500, 400))

  if (!is.null(density_color)) {
    surface_mesh$material <- list(color = density_color)
    shade3d(surface_mesh, meshColor = "vertices", main = "Age", specular = 'black')
  } else {
    shade3d(surface_mesh)
  }

  if (legend) {
    if (is.null(legend_color_sel)) {
      legend_color_sel <- c("grey", "blue", "green", "yellow", "orange", "red", "pink")
    }

    if (is.null(legend_maxi) && is.null(legend_mini)) {
      legend_maxi <- 2100
      legend_mini <- 0
    }

    breaks <- seq(legend_maxi, legend_mini, length.out = length(legend_color_sel))
    legend_labels <- round(breaks, -1)

    legend_colors <- rev(legend_color_sel)

    bgplot3d({
      plot.new()
      legend("topright", title = expression("Density (mg/cm"^3*")"),
             legend = legend_labels,
             fill = legend_colors,
             cex = 1.2,
             bty = "n")

      mtext(title, side = 1, line = 3, cex = 1.4)
    })
  }

  if (hasArg(userMat)) {
    view3d(userMatrix = userMat)
  }
}


#' Plot Cross-Sectional Bone Visualization in 3D
#'
#' Visualizes a 3D cross-section of a bone using surface mesh and internal density
#' (fill) points. Clips the surface mesh at a given axis and value, and overlays a
#' 2D projection of internal density.
#'
#' @param surface_mesh A `mesh3d` object representing the outer surface of the bone.
#' @param surface_colors Optional. A vector of colors for each vertex of the surface mesh. If NULL, uses mesh's own material colors.
#' @param fill_coords A numeric matrix of internal fill point coordinates.
#' @param fill_colors A vector of colors corresponding to fill points.
#' @param slice_axis Character. `'x'`, `'y'`, or `'z'`. Axis along which to slice.
#' @param slice_val Numeric (0 to 1). Relative slice location along selected axis.
#' @param slice_thickness Numeric. Width of the slice (default = 1).
#' @param legend Logical. Optional color bar.
#' @param legend_color_sel Optional character with color gradient
#' @param legend_maxi Numeric. Maximum bone density.
#' @param legend_mini Numeric. Minimum bone density.
#' @param IncludeSurface Logical. Whether to include the clipped surface mesh.
#' @param title Character. Title for the plot.
#' @param userMat Optional. A 4x4 matrix controlling view orientation.
#' @return Generates an `rgl` plot
#' @examples
#' \donttest{
#'   # Download CT scan
#'   url <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.1/test_CT_hip.nii.gz"
#'   scan_filepath <- tempfile(fileext = ".nii.gz")
#'   download.file(url, scan_filepath, mode = "wb")
#'   nifti <- import_scan(scan_filepath)
#'   url2 <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/test_CT_femur.stl"
#'   bone_filepath <- tempfile(fileext = ".stl")
#'   download.file(url2, bone_filepath, mode = "wb")
#'   surface_mesh <- import_mesh(bone_filepath)
#'   landmark_path <- system.file("extdata", "test_femur.mrk.json",
#'                                package = "BoneDensityMapping")
#'   landmarks <- import_lmks(landmark_path)
#'   mapped_coords <- surface_points_template(surface_mesh, landmarks,
#'                                            no_surface_sliders = 100)
#'   mat_peak <- voxel_point_intersect(mapped_coords, nifti)
#'   colored_mesh <- color_mesh(surface_mesh, mapped_coords, mat_peak)
#'   internal_fill <- fill_bone_points(surface_mesh, 3)
#'   internal_density <- voxel_point_intersect(internal_fill, nifti,
#'                                             ct_eqn = "linear",
#'                                             ct_params = c(68.4, 1.106))
#'   internal_colors <- color_mapping(internal_density)
#'   plot_cross_section_bone(colored_mesh, surface_colors = NULL,
#'                           internal_fill, internal_colors, slice_axis = 'x',
#'                           slice_val = 0.5)
#' }
#' @import rgl
#' @import geometry
#' @import concaveman
#' @import sp
#' @importFrom graphics par
#' @export
plot_cross_section_bone <- function(surface_mesh,
                                    surface_colors = NULL,
                                    fill_coords, fill_colors,
                                    slice_axis, slice_val,
                                    slice_thickness = 1,
                                    IncludeSurface = FALSE,
                                    title = "Bone Cross-Section",
                                    userMat = NULL,
                                    legend = TRUE,
                                    legend_color_sel = NULL,
                                    legend_maxi = NULL,
                                    legend_mini = NULL) {
  # check inputs
  stopifnot(nrow(fill_coords) == length(fill_colors))
  if (is.null(slice_axis) || is.null(slice_val)) {
    stop("slice_axis and slice_val must be provided")
  }

  axis_index <- switch(slice_axis, x = 1, y = 2, z = 3,
                       stop("slice_axis must be one of 'x', 'y', or 'z'"))

  # Extract surface coordinates
  surface_coords <- t(surface_mesh$vb[1:3, , drop = FALSE])

  if (is.null(surface_colors) && !is.null(surface_mesh$material$color)) {
    surface_colors <- surface_mesh$material$color
  }

  # Compute slice value from relative position
  coord_val <- min(fill_coords[, axis_index]) +
    slice_val * (max(fill_coords[, axis_index]) - min(fill_coords[, axis_index]))

  open3d()
  par3d(windowRect = c(20, 30, 500, 400))

  shade3d(surface_mesh, color = "gray", alpha = 0.2)

  # Clip surface
  if (IncludeSurface) {
    keep_surface <- surface_coords[, axis_index] <= coord_val
    kept_vertex_indices <- which(keep_surface)

    if (length(kept_vertex_indices) >= 3) {
      clipped_vertices <- surface_mesh$vb[, kept_vertex_indices, drop = FALSE]
      old_faces <- surface_mesh$it

      is_vertex_kept <- rep(FALSE, max(old_faces))
      is_vertex_kept[kept_vertex_indices] <- TRUE

      faces_kept_logical <- colSums(matrix(is_vertex_kept[old_faces], nrow = 3)) == 3
      clipped_faces <- old_faces[, faces_kept_logical, drop = FALSE]

      vertex_map <- integer(max(old_faces))
      vertex_map[kept_vertex_indices] <- seq_along(kept_vertex_indices)
      clipped_faces <- matrix(vertex_map[clipped_faces], nrow = 3)

      clipped_surface_colors <- if (!is.null(surface_colors)) {
        surface_colors[kept_vertex_indices]
      } else {
        "gray"
      }

      clipped_surface_mesh <- surface_mesh
      clipped_surface_mesh$vb <- clipped_vertices
      clipped_surface_mesh$it <- clipped_faces
      clipped_surface_mesh$material <- list()

      clipped_surface_mesh <- Rvcg::vcgClean(clipped_surface_mesh, sel = 0)
      clipped_surface_mesh <- Rvcg::vcgUpdateNormals(clipped_surface_mesh)

      shade3d(clipped_surface_mesh, meshColor = "vertices",
              col = clipped_surface_colors, specular = "black", smooth = TRUE)
    }
  }

  # Slice fill
  keep_fill <- fill_coords[, axis_index] >= (coord_val - slice_thickness / 2) &
    fill_coords[, axis_index] <= (coord_val + slice_thickness / 2)
  fill_slice_coords <- fill_coords[keep_fill, , drop = FALSE]
  fill_slice_colors <- fill_colors[keep_fill]

  if (nrow(fill_slice_coords) >= 3) {
    fill_slice_coords[, axis_index] <- coord_val
    other_axes <- setdiff(1:3, axis_index)
    projected_2d <- fill_slice_coords[, other_axes, drop = FALSE]

    hull_coords <- concaveman(projected_2d)
    tri <- delaunayn(projected_2d)

    keep_tri <- sapply(1:nrow(tri), function(i) {
      idx <- tri[i, ]
      tri_coords <- projected_2d[idx, , drop = FALSE]
      centroid <- colMeans(tri_coords)
      point.in.polygon(centroid[1], centroid[2], hull_coords[, 1], hull_coords[, 2]) > 0
    })

    tri_filtered <- tri[keep_tri, , drop = FALSE]

    if (nrow(tri_filtered) >= 1) {
      mesh <- tmesh3d(
        vertices = t(fill_slice_coords),
        indices = t(tri_filtered),
        homogeneous = FALSE
      )
      mesh$material <- list()
      shade3d(mesh, col = fill_slice_colors, meshColor = "vertices", specular = "black")
    }
  }

  bgplot3d({
    plot.new()
    if (legend) {
      if (is.null(legend_color_sel)) {
        legend_color_sel <- c("grey", "blue", "green", "yellow", "orange", "red", "pink")
      }

      if (is.null(legend_maxi) || is.null(legend_mini)) {
        legend_maxi <- 2100
        legend_mini <- 0
      }
      breaks <- seq(legend_maxi, legend_mini, length.out = length(legend_color_sel))
      legend_labels <- round(breaks, -1)

      legend_colors <- rev(legend_color_sel)

      legend("topright", title = expression("Density (mg/cm"^3*")"),
             legend = legend_labels,
             fill = legend_colors,
             cex = 1.2,
             bty = "n")
    }

    mtext(side = 1, title, line = 3, cex = 1.4)
  })

  if (!is.null(userMat)) {
    view3d(userMatrix = userMat)
  }
}


#' Produce stand alone color bar
#' @param colors String
#' @param mini Numeric
#' @param maxi Numeric
#' @param orientation "horizontal" or "vertical"
#' @param breaks Numeric vector
#' @param title String
#' @param text_size Numeric
#' @param plot Logical
#' @examples
#' colors <- c("darkblue", "blue", "lightblue", "green", "yellow", "red", "pink")
#' color_bar(colors, 0, 2000, breaks = c(0, 500, 1000, 1500, 2000))
#' @importFrom ggplot2 ggplot unit labs guides theme element_text geom_point aes scale_color_gradientn guide_colorbar
#' @importFrom cowplot get_legend
#' @importFrom ggpubr as_ggplot
#' @return ggplot object
#' @export
color_bar <- function(colors, mini, maxi, orientation = "vertical", breaks,
                      title = "", text_size = 11, plot = TRUE) {
  x <- y <- NULL
  # make plot
  z2 <- seq(from = mini, to = maxi, length.out = 100)
  df <- data.frame(x = 1:100, y = 1:100, z2 = z2)
  g <- ggplot(df, aes(x, y)) + geom_point(aes(color = z2))
  g <- g + scale_color_gradientn(colors = colors,
                                 breaks = breaks)
  g <- g + labs(color = title)
  if (orientation == "horizontal") {
    g <- g + guides(color = guide_colorbar(title.position = "top"))
    g <- g + theme(legend.key.size = unit(2.0, "cm"),
                   legend.position = "bottom")
  }
  if (orientation == "vertical") {
    g <- g + theme(legend.key.size = unit(2.0, "cm"),
                   legend.text = element_text(size = text_size))
  }

  legend <- get_legend(g)
  lg <- as_ggplot(legend)
  if (plot == TRUE) {
    lg
  }
  return(lg)
}


#' local significance
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param vertices Matrix
#' @param sig_vals Numeric vector
#' @param changes Numeric vector
#' @param sig_level Numeric. Default 0.05
#' @param dist Numeric. Distance to check for vertices
#' @return Numeric vector
#' @importFrom rdist cdist
rm_local_sig <- function(vertices, sig_vals, changes, sig_level = 0.05,
                         dist) {
  # identify significant values
  sig_inds <- which(sig_vals < sig_level)
  sig_changes <- changes[sig_inds]
  sig_inds_up <- sig_inds[which(sig_changes > 0)]
  sig_inds_down <- sig_inds[which(sig_changes < 0)]

  # check if nearby values are also significant
  sig_vals_updated <- sig_vals
  sig_verts_up <- vertices[sig_inds_up, ]
  sig_verts_down <- vertices[sig_inds_down, ]
  for (i in 1:length(sig_inds_up)) {
    vert <- vertices[sig_inds_up[i], ]
    y <- cdist(vert, sig_verts_up)
    if (length(which(y < dist)) < 2) {sig_vals_updated[sig_inds_up[i]] = 0.1}
  }
  for (i in 1:length(sig_inds_down)) {
    vert <- vertices[sig_inds_down[i], ]
    y <- cdist(vert, sig_verts_down)
    if (length(which(y < dist)) < 2) {sig_vals_updated[sig_inds_down[i]] = 0.1}
  }

  # return vector
  return(sig_vals_updated)
}


#' Sigma beta CT calculations
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param ct_nos Numeric vector. CT numbers from scan
#' @param calibration_type String. Currently only linear supported.
#' @param params Numeric vector. beta and sigma values for calibration eqn
#' @return Vector with estimated density values  in mg/cm^3
ct_calibration <- function(ct_nos, calibration_type, params) {
  # calibration equation
  if (calibration_type == "linear") {
    if (length(params) != 2) {stop("params needs beta and sigma coefficients")}
    density_vals <- (ct_nos - params[1]) / params[2]
  }

  # return
  return(density_vals)
}


#' Reorientate landmarks
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param landmark_path String
#' @param x Integer Value to apply to convert mesh i.e. -1 will mirror x coords
#' @param y Integer Value to apply to convert mesh i.e. -1 will mirror y coords
#' @param z Integer Value to apply to convert mesh i.e. -1 will mirror z coords
#' @return Overwritten landmark file
#' @examples
#' landmark_path <- system.file("extdata", "test_femur.mrk.json",
#'                              package = "BoneDensityMapping")
#' reoriented_landmarks <- reorientate_landmarks(landmark_path)
#' @importFrom utils read.table write.table
#' @importFrom jsonlite read_json write_json
#' @export
reorientate_landmarks <- function(landmark_path, x = 1, y = 1, z = 1) {
  # Check file extension
  file_ext <- tools::file_ext(landmark_path)

  if (file_ext == "fcsv") {
    # Handle FCSV format
    header <- readLines(landmark_path, n = 3)
    lmks <- read.table(landmark_path, sep = ",", skip = 3, header = TRUE)

    # Apply mirroring
    lmks[, "x"] <- lmks[, "x"] * x
    lmks[, "y"] <- lmks[, "y"] * y
    lmks[, "z"] <- lmks[, "z"] * z

    # Overwrite file with header and new landmark coordinates
    writeLines(header, con = landmark_path)
    write.table(lmks, file = landmark_path, append = TRUE, sep = ",",
                row.names = FALSE, col.names = TRUE)

  } else if (file_ext == "json") {
    # Handle JSON format
    data <- read_json(landmark_path)  # no simplifyVector!

    for (i in seq_along(data$markups[[1]]$controlPoints)) {
      coords <- as.numeric(data$markups[[1]]$controlPoints[[i]]$position)
      data$markups[[1]]$controlPoints[[i]]$position <- c(coords[1] * x,
                                                         coords[2] * y,
                                                         coords[3] * z)
    }

    write_json(data, landmark_path, auto_unbox = TRUE, pretty = TRUE)

  } else {
    stop("Unsupported file type. Only .fcsv and .json formats are supported.")
  }
}
