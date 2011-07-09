/*---------------------------------------------------------------------------##
##  Library:
##      galosh::profillic
##  File:
##      Profillic.hpp
##  Author:
##      D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
##  Description:
##      Global #defines for the Profillic software for estimating parameters of profile HMMs.
##
#******************************************************************************
#*
#*    This file is part of profillic, a suite of programs for estimating parameters of
#*    Profile HMMs.  Please see the document CITING, which should have been
#*    included with this file.  You may use at will, subject to the license
#*    (Apache v2.0), but *please cite the relevant papers* in your documentation
#*    and publications associated with uses of this library.  Thank you!
#*
#*    Copyright (C) 2008, 2011 by Paul T. Edlefsen, Fred Hutchinson Cancer
#*    Research Center.
#*
 *    Licensed under the Apache License, Version 2.0 (the "License");
 *    you may not use this file except in compliance with the License.
 *    You may obtain a copy of the License at
 *    
 *        http://www.apache.org/licenses/LICENSE-2.0
 *    
 *    Unless required by applicable law or agreed to in writing, software
 *    distributed under the License is distributed on an "AS IS" BASIS,
 *    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *    See the License for the specific language governing permissions and
 *    limitations under the License.
#*****************************************************************************/

#if     _MSC_VER > 1000
#pragma once
#endif

#ifndef __GALOSH_PROFILLIC_HPP__
#define __GALOSH_PROFILLIC_HPP__

// ProfileTrainer
#define GALOSH_INLINE_TRAIN
#define GALOSH_INLINE_TRAIN_PERCENT_CHANGE inline

// ProfileGibbs
#define GALOSH_INLINE_GIBBS

// SimulationStudy
#define GALOSH_INLINE_SIMULATIONSTUDY_START

#endif // __GALOSH_PROFILLIC_HPP__
