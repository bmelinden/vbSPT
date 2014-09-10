 vbSPT is a program suite for analysis of single particle diffusion
 trajectories, where the diffusion constants switch randomly according
 to a discrete Markov process. The program runs on Matlab, but uses
 compiled C-code to speed up the most computer intensive loops.  The
 latest version of the software as well as a forum for discussion and
 questions can be found at 'sourceforge.net/projects/vbspt/'.

 ========================================================================= 
 Copyright (C) 2013 Martin Lindén, Fredrik Persson, and Johan Elf
 
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

 This product includes software developed by Martin Lindén (in HMMcore/), 
 and by Jan-Willem van de Meent (in external/), see copyright in individual 
 files.
 =========================================================================
v 1.1.1 (2014-10-11): minor bugfix release
 - corrected a sign error in the KL divergence of the initial state vector 
  (the faulty terms was gammaln of the total prior weight, and so very 
   unlikely to make any practical difference).
 - corrected some internal function names to be consistent with file names.
 - made VB3_readData handle input files without the .m ending
--------------------------------------------------------------------------
v 1.1 (2013-11-12): code speed-up
 - better representation of data leads to a significantly faster
   execution (~10-fold on our test data).
 - A new parameterization of the transition rates, allowing more
   flexible priors for advanced users (See manual). For non-advanced
   users, the old priors (Nat. Meth 2013 paper) are still default.
 - updated VB3_displayHMMmodel and VB3_getTrjStats to the new model format.
 - deleted VB3_reorder.m (obsolete, replaced by VB3_sortModel.m)
--------------------------------------------------------------------------
 v 1.0.1 (2013-05-31): bugfix release
 - getDwellTRJ.m : consistent function name in m-file, and better handling 
   trajectories with not all states present.
 - VB3_HMManalysis.m : added functionaly to create non-existent target
   folders.
 - VBviterbi_log.c   : fixed bug in viterbi algorithm.
--------------------------------------------------------------------------
