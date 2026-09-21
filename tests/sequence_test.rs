use frugal::sequence::{amino_letter, amino_num, imin, max_fr, rframe};
use std::os::raw::{c_char, c_int};

#[test]
fn test_imin() {
    unsafe {
        assert_eq!(imin(3, 5), 3);
        assert_eq!(imin(5, 3), 3);
        assert_eq!(imin(-1, 0), -1);
    }
}

#[test]
fn test_amino_num_letter_roundtrip() {
    unsafe {
        let letters = b"ACDEFGHIKLMNPQRSTVWY";
        for (i, &c) in letters.iter().enumerate() {
            assert_eq!(amino_num(c as c_char), i as c_int);
            assert_eq!(amino_letter(i as c_int) as u8, c);
        }
        assert_eq!(amino_num(b'Z' as c_char), -1);
        assert_eq!(amino_letter(20) as u8, b'X');
    }
}

#[test]
fn test_rframe() {
    unsafe {
        assert_eq!(rframe(0, 10), 3);
        assert_eq!(rframe(1, 10), 2);
    }
}

#[test]
fn test_max_fr() {
    unsafe {
        assert_eq!(max_fr(3, 2, 1), 0);
        assert_eq!(max_fr(1, 3, 2), 1);
        assert_eq!(max_fr(1, 2, 3), 2);
        // ties go to 2
        assert_eq!(max_fr(1, 1, 1), 2);
    }
}
