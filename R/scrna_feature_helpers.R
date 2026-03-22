#' CSS vascular gate
#'
#' @param R Vascular proportion.
#' @param tau Threshold.
#'
#' @return Numeric vector.
#' @noRd
css_gate_fun <- function(R, tau) {
  pmin(1, R / tau)
}

#' Sender prior
#'
#' @param sender Sender cell type label.
#'
#' @return Numeric scalar.
#' @noRd
pi_sender <- function(sender) {
  s <- tolower(as.character(sender))
  if (s %in% c("vascular")) return(1.00)
  if (grepl("immune|microglia|macrophage", s)) return(0.80)
  if (grepl("stem|progenitor|npc", s)) return(0.70)
  if (grepl("glia|astro|oligo", s)) return(0.30)
  if (grepl("neuron", s)) return(0.20)
  0.50
}

#' Split complex receptor names
#'
#' @param x Character scalar.
#'
#' @return Character vector.
#' @noRd
split_complex <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  parts <- unlist(strsplit(as.character(x), "_", fixed = TRUE))
  parts <- parts[parts != ""]
  toupper(parts)
}

#' Standardize CSS table to group-level vascular ratios
#'
#' @param dt CSS table.
#'
#' @return A data.table with columns `groupAB` and `R_v`.
#' @noRd
standardize_css <- function(dt) {
  dt <- data.table::as.data.table(dt)
  cols <- names(dt)
  req <- c("groupAB", "cell_type_major", "mean_pct")

  if (!all(req %in% cols)) {
    stop(
      "[CSS] Expect columns: ", paste(req, collapse = ", "),
      "\n  Available columns: ", paste(cols, collapse = ", "),
      call. = FALSE
    )
  }

  dt2 <- data.table::copy(dt)
  dt2[, cell_type_major := tolower(as.character(cell_type_major))]

  v <- dt2[
    grepl("vasc", cell_type_major),
    list(R_v = mean(mean_pct, na.rm = TRUE)),
    by = "groupAB"
  ]

  if (nrow(v) == 0L) {
    stop(
      "[CSS] No vascular row found in cell_type_major. Unique cell_type_major = ",
      paste(unique(dt2$cell_type_major), collapse = ", "),
      call. = FALSE
    )
  }

  if (max(v$R_v, na.rm = TRUE) > 1.5) {
    v[, R_v := R_v / 100]
  }

  v[]
}

#' Prepare CIS table
#'
#' Keeps receiver = vascular when present, adds CSS gate, sender prior,
#' and computes the normalized interaction score used in v1.0.
#'
#' @param dt CIS table.
#' @param css_grp Output of `standardize_css()` plus `Cp`.
#' @param eps Small constant for numerical stability.
#'
#' @return A prepared data.table with `S_int`.
#' @noRd
prep_cis <- function(dt, css_grp, eps = 1e-12) {
  dt <- data.table::as.data.table(dt)
  css_grp <- data.table::as.data.table(css_grp)

  cols <- names(dt)
  req <- c("groupAB", "sender", "n_interactions", "sum_weight")

  if (!all(req %in% cols)) {
    stop(
      "[CIS] Expect columns: ", paste(req, collapse = ", "),
      "\n  Available columns: ", paste(cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (!all(c("groupAB", "Cp") %in% names(css_grp))) {
    stop("[CIS] css_grp must contain columns: groupAB, Cp", call. = FALSE)
  }

  d <- data.table::copy(dt)

  if ("receiver" %in% cols) {
    d[, receiver := tolower(as.character(receiver))]
    d <- d[receiver == "vascular"]
  }

  d[, sender := tolower(as.character(sender))]

  d <- merge(
    d,
    css_grp[, c("groupAB", "Cp"), with = FALSE],
    by = "groupAB",
    all.x = TRUE
  )

  d[is.na(Cp), Cp := 0]

  d[, pi_s := vapply(sender, pi_sender, numeric(1))]

  d[, n_tilde := (n_interactions + eps) /
      (sum(n_interactions + eps, na.rm = TRUE) + eps),
    by = "groupAB"
  ]
  d[, w_tilde := (sum_weight + eps) /
      (sum(sum_weight + eps, na.rm = TRUE) + eps),
    by = "groupAB"
  ]

  d[!is.finite(n_tilde), n_tilde := 0]
  d[!is.finite(w_tilde), w_tilde := 0]

  d[, I_base := 0.5 * (n_tilde + w_tilde)]
  d[, c("n_tilde", "w_tilde") := NULL]

  d[, S_int := Cp * I_base * pi_s]

  d[]
}

#' Prepare LRS table
#'
#' Keeps target = vascular when present, uppercases ligand/receptor symbols,
#' and computes base edge weight `A0 = weight^gamma`.
#'
#' @param dt LRS table.
#' @param gamma Exponent applied to raw edge weight.
#'
#' @return A prepared data.table.
#' @noRd
prep_lrs <- function(dt, gamma = 1) {
  dt <- data.table::as.data.table(dt)
  cols <- names(dt)
  req <- c("groupAB", "source", "ligand", "receptor", "weight")

  if (!all(req %in% cols)) {
    stop(
      "[LRS] Expect columns: ", paste(req, collapse = ", "),
      "\n  Available columns: ", paste(cols, collapse = ", "),
      call. = FALSE
    )
  }

  d <- data.table::copy(dt)

  d[, source := tolower(as.character(source))]

  if ("target" %in% cols) {
    d[, target := tolower(as.character(target))]
    d <- d[target == "vascular"]
  }

  d[, ligand := toupper(as.character(ligand))]
  d[, receptor := toupper(as.character(receptor))]
  d[, weight := as.numeric(weight)]
  d[is.na(weight), weight := 0]

  d[, A0 := weight^gamma]

  d[]
}
