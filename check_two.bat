SET thresh=20
SET err=0.03
SET runflags=-t -v -i

cargo build --release --features "valimaki2"
move target\release\rust_overlaps.exe rust_overlaps_v.exe
cargo build --release --features "kucherov"
move target\release\rust_overlaps.exe rust_overlaps_k.exe
rust_overlaps_v.exe TEST.fasta TEST_v.txt %err% %thresh% %runflags%
rust_overlaps_k.exe TEST.fasta TEST_k.txt %err% %thresh% %runflags%
powershell -command "& {&'diff' (cat TEST_v.txt) (cat TEST_k.txt) > TEST_diff.txt}"
TIMEOUT 10