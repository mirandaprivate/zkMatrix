/// Two-tier commitment to a sparse matrix a.
/// 
/// * The first tier is to commit to each row(column) of the matrix.
/// * The second tier is to commit to the first tier commitments
///     by using the AFGHO commitment.
/// 
/// This is equivalent to:
/// $$
///     C_a = \sum_{i=1}^{n} a_{ij} e(G_i, H_j).
/// $$
/// where
/// $$
///     a_{ij} \in \mathbb{Z}_p,
///     G_i \in G_1, 
///     H_j \in G_2.
/// $$
/// 
/// Note that this commitment can be computed parallelly.
/// We use multi-threading to compute the commitment.
/// 