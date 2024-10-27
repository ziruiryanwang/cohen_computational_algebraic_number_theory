pub trait GroupElement: Clone + PartialEq {
    fn identity() -> Self;
    fn mul(&self, other: &Self) -> Self;
    fn inverse(&self) -> Self;
}
