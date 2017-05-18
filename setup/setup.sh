#
# setup.sh
#
# setup script for moriiism/mchanlc
#
# 2017.05.18 M. Morii
#

#
# This script must be loaded by "source" command of bash,
# before using moriiism/mchanlc
#

##
## mchanlc
##

export MITOOL=/home/morii/work/github/moriiism/mitool
export LD_LIBRARY_PATH=/soft/root/6.08.02/lib:/soft/gsl/2.3/lib:${LD_LIBRARY_PATH}
export PATH=/soft/root/6.08.02/bin:${PATH}

alias root="root -l"

#
# HEADAS
#

HEADAS_VER=6.16
export HEADAS=/soft/heasoft/heasoft-${HEADAS_VER}/x86_64-unknown-linux-gnu-libc2.12
source $HEADAS/headas-init.sh

################################################################
################################################################
############ Don't Edit Below. #################################

##
## change the terminal title
##

termtitle="mchanlc"
PROMPT_COMMAND='echo -ne "\033]0;${termtitle}\007"'
