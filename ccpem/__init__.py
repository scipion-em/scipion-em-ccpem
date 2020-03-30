# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

"""
This package contains the protocols for
manipulation of atomic struct objects
"""

import os
import pwem
import pyworkflow.utils as pwutils
from .bibtex import _bibtexStr

_logo = 'ccpem.png'

class Plugin(pwem.Plugin):
    _homeVar = 'CCPEM_HOME'

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar('CCPEM_HOME', 'ccpem-1.4.1')

    @classmethod
    def defineBinaries(cls, env):
        ccpem_commands = [('./install_ccpem.sh', ['setup_ccpem.sh'])]

        env.addPackage('ccpem', version='1.4.1',
                       url='https://www.ccpem.ac.uk/downloads/ccpem_distributions/ccpem-1.4.1-linux-x86_64.tar.gz',
                       commands=ccpem_commands,
                       default=True)

    @classmethod
    def getEnviron(cls, ccpemFirst=True):
        """ Create the needed environment for Xmipp programs. """
        environ = pwutils.Environ(os.environ)
        pos = pwutils.Environ.BEGIN if ccpemFirst else pwutils.Environ.END
        environ.update({
            'PATH': cls.getHome('bin'),
            'CLIB': cls.getHome('lib'),
            'CLIBD': cls.getHome('lib/data'),
            'CLIBDMON': cls.getHome('lib/data/monomers'),
        }, position=pos)

        return environ


