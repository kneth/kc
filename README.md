# kc - A Kinetic Compiler

Source code, test files, and documentation in this repository go back
to mid-1990s. With only a few modifications, the files are left as
they were. The authors' addresses have been removed since they have
changed. Moreover, the license has changed from GNU General Public
License v2 to GNU General Public License v3. Furthermore this preamble
has been added.

Recently (November 2020) the source code has been compiled on MacOS
v10.15 using GCC.

The publication of kc has been a wish for some time. We hope it can be
useful to some, and asumement for others.

## Introduction

This is the README file for the kc project. The kc program is a
"kinetic compiler". This means it is able to transform chemical
equations into simulation programs.

The motivation of the compiler and the use is documented in one of the
authors' master thesis (_Chemical waves in reaction-diffusion systems:
a numerical study_, K. Geisshirt, University of Copenhagen, May 1994).

The newest edition of the user's manual can be founded in
`docs/kc-man.tex`, and the programmer's manual in `docs/kc.tex`. For
users who do not have LaTeX, please write to the author to obtain
either a hard-copy or a postscript version.

In addition to the user's manual there exists a directory with
examples of input to kc. The directory is named `test` (some of our test
examples are found here).

No program is bug free, and the known bugs in kc are documented in the
file `KnownBugs.md`. Furthermore, the program is not complete (which
program is?), and a list of future enchancement is found in the file
`TODO.md`.


## The history of kc
*  April 1995: Version 1.05 released. A number of bugs were 
   discovered in version 1.0 and they have been fixed. Some
   cleaning up has been done.

   This version was developed at KIKU and RUC.

*  December 1994: Version 1.00 released. The compiler is now running
   very fine. It works very well together with a backend which is 
   solving ordinary differential equations. The following code 
   generators are supported in this version:
     1. Waves (dymanical simulations of reaction-diffusion systems).
     2. Eigen (computing various properties of a given system of ODEs).
     3. KNcont (continuation program written in Pascal).
     4. KCI (numerical solution of ODEs/dynamical simulation of 
        homogenous chemical reactions).

*  July 1994: Version 0.99 released (only for internal use). A number
   of bugs have been fixed. A solver of ordinary differential
   equations called KGode has been developed. The solver is written in
   ANSI C and have been ported to at least to four platforms: HP-UX,
   ConvexOS, Linux, and MS-DOS. There exists an automatic steplength
   controller. The solver implements three numerical schemes:
     1. 4th order Runge-Kutta [1]. 
     2. Calahan's method [2].
     3. Generalised Runge-Kutta [1, 3].

   Various modes have been developed and old ones have been extended.
   But in total, the following modes are supported now:
     1. Kin (dynamical simulations - requires a Pascal compiler).
     2. Waves (dynamical simulations of reaction-diffusion systems,
        simulation programs written by K. Geisshirt).
     3. KNcont (continuations - require a Pascal compiler).
     4. Eigen (calculation of eigenvalues and -vectors of the Jacobian
        matrix). 
     5. KGode (dynamical simulations).
     6. Stoc (dynamical simulations using a stochastic approach [5]).
     7. IScont (also called CONT - a continuation package in
        Fortran-77). 

   The grammar has been extended a little so now it is possible of
   read and store string constants. The grammar for expressions has
   been modified a bit in order to fix a major bug. 

   This version has been developed at KIKU.

*  December 1993: Version 0.50 released. It mainly fixes of a number of bugs.

*  September 1993: Version 0.25 released. Only the following modes are
   supported in this version:
     1. Kin
     2. Waves
     3. KNcont
     4. Eigen
   Other modes in previous versions were never used, and therefore the
   author stopped supporting them.

   A lot of cleaning up has been done. K. Nielsen has done a lot of
   testing, which has discovered many minor bugs and some major ones.
   These errors had been corrected. 

*  May 1993: Version 0.20 released. Autocatalytic reactions can be
   used, and the code generators are supporting parameters. Bugs in
   version 0.10 have been corrected. This version is developed at KIKU
   in Copenhagen. The grammar is changed, so the following features
   are supported: 
     1. parameters

   The following code generators have been added:
     1. Waves (A PDE solver by O. Jensen et al.).
     2. Keld Nielsen's continuation program (KNcont).
     3. A mode which calculates the Jacoby matrix, its eigenvalues and
   eigenvectors. 

   The code generator to CONT is now full operative, i.e. it supports
   the use of parameters. 

*  October 1992: Version 0.10 released. This version supports KIN,
   CONT, and Dalimil Snita's Chemical Meta Language. This version is
   developed at KIKU in Copenhagen. The grammar is now supporting 
     1. ordinary differential equations
     2. constraints 

   The program is now not using so much memory as previous version.

*  September 1992: Version 0.00 released. This version only supports
   KIN and CONT. This version is the initial guess of the system. It
   was developed at VSCHT in Prague.


## Installation and running

The installation is very easy! There exists a small script called
kc-inst which does most of the job. The script must be supplied with
two arguments, namely a directory name and a platform. The directory
name is a prefix to the directories which you want kc to "live". The
platform is at the moment one of the following: 

* Generic GNU C-compiler
* Linux 
* DECstations running Ultrix 
* Silicon Graphics 
* Convex running ConvexOS 
* HP-UX v7 
* IBM RS/6000 running AIX v3.2

Run the script with no arguments and you will get some help. The
script will generate a script called kci, which is a front-end to the
kci-mode. 

An example: If you want to install the system on a Linux machine 
in the directories /usr/local you should type:

```sh
kc-inst /usr/local LINUX
```

In order to use the DOS version, you have to have an Intel 80386 (or
higher) computer. You also have to install the djgpp package, and this
package is surposed to be used during all your simulations. The
makefile you should used is `src/Makefile.DOS`, and the file
`Scripts/kci.bat` is a front-end to the `kci-mode`.

There can be a problem with newer Unices, because they have support
for your natioonal language. Therefore, check the enviroment variable
called `LANG` to see if it is set to `C`.

Another thing is that you should create a new directories before
installing. Let us assume you want to use the directory `DIR` as prefix
(as mentioned above). You should them create the directories `DIR/num`
and `DIR/bin`.

## Legal issues

The program is copyrighted in the sense of GNU General Public License
v3. The authors cannot be responsible for any looses the programs in
this package may produce.

If the program is used to published scientific results, authors are
asked to make a footnote or better a reference to kc. A the moment the
best reference is [4], but this may change in the future.

The package comes under the GNU General Public License v3 (GPLv3), and GPLv3
is supplied in form of the file `LICENSE`. Please read this file.


## Support

The program is NOT supported by the authors, but you are welcome to
create an issue.


## Acknowledgements

The authors wish to thank the following people:

* M. Marek - ideas to the first version.
* P.G. Sorensen - do. and advices.
* F. Hynne - same as P.G. Sorensen.
* and the first users (J. Wang, A. Nagy, etc.).

## Abbreviations

* VSCHT:  Prague Institute for Chemical Engineering.
* KIKU: Institute of Chemistry, University of Copenhagen.
* KIN: KINetic compiler and simulator (by P.G. Sorensen).
* CONT: CONTinuation program (by I. Screiber).
* EIGEN: Mode calculating the jacoby matrix and its eigenvectors and eigenvalues.
* KNcont: Keld Nielsen's continuation programs.
* KGode: An ODE solver written by Keld Nielsen and Kenneth Geisshirt.
* Stoc: A dynamical simulator using a stochastic scheme.
* GPL: GNU General Public License.
* RUC: Roskilde University
