#!/usr/bin/env python

#    xinterval.py
#    Copyright (C) 2006, Simon David Burton and The Australian Natonal University.
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


import sys

if 'optimize' in sys.argv:
    sys.argv.remove('optimize')
    from optimize import Interval, get_default_prec, set_default_prec
else:
    from mpfi import Interval, get_default_prec, set_default_prec, DisorderException




