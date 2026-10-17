#!/bin/bash
# Download the publisher's Supplementary Information, Reporting Summary, Peer Review File, Supplementary Tables 1-7 and the Source Data workbooks
# of Gaborieau et al. 2024 (13 files, 8.7 MB) into docs/paper/supplement/. They are not redistributed in the repository.
#   scripts/download_supplement.sh
cd "$(dirname "$0")/.." || exit 1
D=docs/paper/supplement; mkdir -p $D
B="https://static-content.springer.com/esm/art%3A10.1038%2Fs41564-024-01832-5/MediaObjects/41564_2024_1832_MOESM"
# 1 Supplementary Information, 2 Reporting Summary, 3 Peer Review File, 4 Supplementary Tables 1-7, 5-13 Source Data (Fig 1, 2, 3, 4, 5, 6, ED6, ED7, ED9)
for i in 1 2 3; do [ -s $D/MOESM$i.pdf ] || curl -fsS -o $D/MOESM$i.pdf "${B}${i}_ESM.pdf"; done
for i in 4 5 6 7 8 9 10 11 12 13; do [ -s $D/MOESM$i.xlsx ] || curl -fsS -o $D/MOESM$i.xlsx "${B}${i}_ESM.xlsx"; done
ls -l $D | grep MOESM
