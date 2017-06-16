use fmt;

#[inline]
pub fn for_reversed_string(id : usize) -> bool {
    id % 2 == 1
}

/*
only needed if reversals are enabled.
For each input string, two unique text entries are appended, one forwards and one backwards.
The orientation of these strings corresponds with the parity of their ID's
(an internal representation of the strings) with NORMAL strings having EVEN IDs and REVERSED strings having ODD parity.
Note that this is "normal" and "reversed" considering the inverse reading direction of the index.
ie: given a string "XYZ", this will append the "normal" (ZYX) to the text with ID 0 and (XYZ) to the text with ID 1.

This function simply finds the ID of the strings with the given ID such that the two
strings are from the same input string, but lie in opposite directions
*/
#[inline]
pub fn companion_id(id : usize, reversals : bool) -> usize{
    assert!(reversals);
    if for_reversed_string(id) {id-1} else {id+1}
}

pub fn relative_orientation(id_a : usize, id_b : usize, reversals : bool) -> Orientation {
    if reversals && (id_a + id_b)%2 == 1{
        Orientation::Reversed
    } else {
        Orientation::Normal
    }
}

// Normal refers to a solution where both strings are NOT reversed
// Reversed refers to a solution where A is normal and B is reversed
#[derive(Hash,PartialEq, Eq, Debug, PartialOrd, Ord, Clone)]
pub enum Orientation{
    Normal,
    Reversed,
}

impl fmt::Display for Orientation {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let s = if *self == Orientation::Normal{"N"} else {"I"};
        write!(f, "{}", s)
    }
}

