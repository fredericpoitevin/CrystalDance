#!/bin/bash
BIN=../source
# PHONON
$BIN/phonon t 0.0 0.0 0.0  > o_pt
$BIN/phonon 0 0.0 0.0 0.01 > o_p0
$BIN/phonon 1 0.0 0.0 0.1  > o_p1
$BIN/phonon 2 0.0 0.0 0.2  > o_p2
$BIN/phonon 3 0.0 0.0 0.3  > o_p3
$BIN/phonon 4 0.0 0.0 0.4  > o_p4
$BIN/phonon 5 0.0 0.0 0.5  > o_p5
$BIN/phonon 6 0.3 0.3 0.3  > o_p6
#
# INCOH
$BIN/incoh 6 13 2 > o_in
#
# COH
$BIN/coh 6 5 x 2 > o_ch
#
# BEAD
$BIN/bead 1 > o_bd

