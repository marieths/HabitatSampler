# HaSa - HabitatSampler
#
# Copyright (C) 2020  Carsten Neumann (GFZ Potsdam, carsten.neumann@gfz-potsdam.de)
#
# This software was developed within the context of the project
# NaTec - KRH (www.heather-conservation-technology.com) funded
# by the German Federal Ministry of Education and Research BMBF
# (grant number: 01 LC 1602A).
# The BMBF supports this project as research for sustainable development (FONA); www.fona.de.
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

################################################################################
clip <- function(raster, shape) {
  raster::rasterOptions(progress = "text")
  a1_crop <- raster::crop(raster, shape)
  step1 <- raster::rasterize(shape, a1_crop)
  step1 <- raster::reclassify(step1, c(1, 200, 1))
  a1_crop * step1
}
