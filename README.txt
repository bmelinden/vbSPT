 vbSPT is a program suite for analysis of single particle diffusion
 trajectories, where the diffusion constants switch randomly according
 to a discrete Markov process. The program runs on Matlab, but uses
 compiled C-code to speed up the most computer intensive loops.  The
 latest version of the software as well as a forum for discussion and
 questions can be found at 'sourceforge.net/projects/vbspt/'.

 ========================================================================= 
 Copyright (C) 2012 Martin Lindén, Fredrik Persson, and Johan Elf
 
 E-mail: 
 bmelinden@gmail.com, freddie.persson@gmail.com, johan.elf@gmail.com
 =========================================================================

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or any
 later version.  This program is distributed in the hope that it will
 be useful, but WITHOUT ANY WARRANTY; without even the implied
 warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
 the GNU General Public License for more details.
 
 If you use vbSPT please cite: "Person F,Lindén M, Unoson C, Elf J,
 Extracting intracellular reaction rates from single molecule tracking
 data, Nature Methods 10, 265–269 (2013). doi:10.1038/nmeth.2367"
 
 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.

 This product includes software developed and copyrighted by Martin
 Lindén (in HMMcore/, see copyright in individual files).

 =========================================================================
 v 1.1 -- beta 1

 Changes since v 1.0 (Nat. Meth. 2013)
 - bugfix in HMMcore/getDwellTRJ.m
 - reparameterization of transition rates, to make the choice of
   priors more transparent (see documentation).
 - bugfix in Viterbi algorithm

 
