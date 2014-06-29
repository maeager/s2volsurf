# README.TXT
#
# Copyright 2012-2014 David G. Barnes
#
# This file is part of S2VOLSURF.
#
# S2VOLSURF is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# S2VOLSURF is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with S2VOLSURF.  If not, see <http://www.gnu.org/licenses/>. 
#
# We would appreciate it if research outcomes using S2VOLSURF would
# provide the following acknowledgement:
#
# "Three-dimensional visualisation was conducted with the S2PLOT
# progamming library"
#
# and a reference to
#
# D.G.Barnes, C.J.Fluke, P.D.Bourke & O.T.Parry, 2006, Publications
# of the Astronomical Society of Australia, 23(2), 82-93.
#
# $Id: README.TXT 78 2012-11-26 22:41:06Z barnesd $
#
#

S2VOLSURF contains utility source code and command-line programs for:

o creating compressed, 8-bit-per-pixel volumes ("xrw volumes") from
  NIFTI images and stacks of TGA images (nifti2xrw, tgastack2xrw);

o displaying basic information about xrw volumes including a
  histogram of voxel values (xrwinfo, xrwhisto); 

o converting xrw volumes to mosaiced PNG files and to simple
  ASCII-format point sets (xrw2pngmos, xrw2points);

o displaying one or more STL-format surface files (s2stl.c - 
  source code must be modified to choose your own STL file/s); 

o extracting an STL-format isosurface from an XRW-format volume
  (xrw2stl.c); 

o displaying the coordinate range of one or more OBJ-format
  surfaces (objrange.c); and

o generating interactive 3-d renderings of xrw volumes, overlaid
  with OBJ-format surfaces, and creating PRC and PDF versions of the
  visualisation (xrw2pdf.c).

See libxrw.c to understand the xrw volume format and find functions
for loading or creating, normalising, and saving xrw volumes.

See libobj.c for functions related to loading and rendering OBJ-format
surfaces.

S2PLOT (version 3.2.1 or higher, http://astronomy.swin.edu.au/s2plot) 
and zlib are required.

PGPLOT (http://www.astro.caltech.edu/~tjp/pgplot) is required to build
xrwhisto.

NIFTICLIB (http://niftilib.sourceforge.net) is required to build
nifti2xrw.

The standard PNG library is required to build xrw2pngmos.

Build the tools using the "build.csh" script.  It may need editing for
appropriate library and include directories for your system. The build
script calls the S2PLOT build scripts, so if S2INSTALLPATH is set,
binaries will be moved to that location after building.

Most tools will print usage information if executed with no arguments,
or invalid arguments.

See EXAMPLE.TXT for an example run-through of xrw2pdf to create a 3-d 
PDF figure (volume rendering + isosurface) from a NIFTI image.

*  Please note that other than the relatively complete and functional
*  xrw2pdf, the tools in the S2VOLSURF kit are best considered
*  "starting points" and may need additional coding to produce the
*  exact output required.  For example, s2stl.c must be modified to
*  refer to your own input STL file/s for it to work.


