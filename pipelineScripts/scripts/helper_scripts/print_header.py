#!/usr/bin/env python3
# This script prints any comment in a structured and pretty way.
# Example usage: announcement () { print_header.py "$1" "H1"; }
# announcement "Step 1"

from __future__ import print_function
import sys

comm = sys.argv[1]
heading = sys.argv[2]

max_len=90

if heading == "H1":
  delim = "#"
  side = "#"
  fill = " "

if heading == "H2":
  delim = "-"
  side = "|"
  fill = " "

if heading == "H3":
  delim = " "
  side = " "
  fill = "-"

if heading == "#":
  delim = " "
  side = " "
  fill = " "

if heading != "#":
  print('\n\n' + delim*120)
  print(side + " "*118 + side)

if len(comm) > max_len:
  cut=comm.split(" ")
  line = ""
  for word in cut:
    if (len(line) + 1 + len(word)) > max_len:
      edge1 = round((120 - len(line))/2) - 2
      edge2 = 120 - edge1 - len(line) - 4
      print(side + fill*edge1 + " " + line + " " + fill*edge2 + side)
      line = word
    else:
      line = line + " " + word
  edge1 = round((120 - len(line))/2) - 2
  edge2 = 120 - edge1 - len(line) - 4
  print(side + fill*edge1 + " " + line + " " + fill*edge2 + side)
else:
  edge1 = round((120 - len(comm))/2) - 2
  edge2 = 120 - edge1 - len(comm) - 4
  print(side + fill*edge1 + " " + comm + " " + fill*edge2 + side)

if heading != "#":
  print(side + " "*118 + side)
  print(delim*120 + '\n\n')


