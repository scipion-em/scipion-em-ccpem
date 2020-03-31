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
import os
import sys

from pwem.protocols import EMProtocol
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import (EnumParam, PointerParam, BooleanParam, IntParam)
from pyworkflow.utils.path import moveFile
from pwem.emlib.image import ImageHandler
from pwem.objects import Volume
from ccpem import Plugin

class ProtFDR(EMProtocol):
    """Applies segmentation by False Discovery Rate"""
    _label = 'fdr'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVol', PointerParam, pointerClass="Volume",
                       label='Volume to segment:', allowsNull=False)
        form.addParam('noiseBox', IntParam, default=20,
                       label='Noise (px):',
                       help='This boxsize should fit in the space between the protein and the volume box. '
                            'After running the protocol, make sure the box is correct from extra/diag_image.pdf')
        form.addParam('useHalfMaps', BooleanParam, default=False,
                       label='Use half maps:')
        form.addParam('method', EnumParam, default=0, choices=['FDR Benjamini-Yekutieli','FDR Bejamini-Hochberg',
                                                               'FWER Holm', 'FWER Hochberg'],
                       label='Method:')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('fdrStep')

    def fdrStep(self):
        args=Plugin.getHome('lib/py2/FDRcontrol.pyc')

        ih = ImageHandler()
        if self.useHalfMaps.get():
            halfmaps = self.inputVol.get().getHalfMaps().split(',')
            fn1=halfmaps[0]
            fn2=halfmaps[1]
            fn1mrc = self._getTmpPath('vol1.mrc')
            fn2mrc = self._getTmpPath('vol2.mrc')
            ih.convert(fn1,fn1mrc)
            ih.convert(fn2,fn2mrc)
            args+=' --em_map tmp/vol1.mrc --halfmap2 tmp/vol2.mrc'
        else:
            fnVol = self._getTmpPath('vol.mrc')
            ih.convert(self.inputVol.get(),fnVol)
            args+=' --em_map tmp/vol.mrc'
        args+=' --testProc rightSided'
        args+=' --window_size %d'%self.noiseBox.get()
        if self.method.get()==0:
            args+=' -method BY'
        elif self.method.get()==1:
            args+=' -method BH'
        elif self.method.get()==2:
            args+=' -method Holm'
        elif self.method.get()==3:
            args+=' -method Hochberg'

        self.runJob('ccpem-python',args, cwd=self._getPath())
        moveFile(self._getPath('diag_image.pdf'),self._getExtraPath('diag_image.pdf'))
        if self.useHalfMaps.get():
            moveFile(self._getPath('vol1_confidenceMap.mrc'), self._getExtraPath('vol_confidenceMap.mrc'))
        else:
            moveFile(self._getPath('vol_confidenceMap.mrc'),self._getExtraPath('vol_confidenceMap.mrc'))

        fhSummary = open(self._getPath('summary.txt'),'w')
        fhLog = open(self._getPath('logs/run.stdout'))
        for line in fhLog.readlines():
            if line.startswith('Calculated map threshold'):
                fhSummary.write(line)
        fhLog.close()
        fhSummary.close()

    # --------------------------- UTILS functions ------------------
    def _validate(self):
        errors = []
        if self.useHalfMaps.get():
            if not self.inputVol.hasValue():
                errors.append("You must provide an input volume")
                return errors
            halfmaps = self.inputVol.get().getHalfMaps().split(',')
            if len(halfmaps)!=2:
                errors.append("The input volume must have two associated half maps.")
        return errors

    def _summary(self):
        summary=[]
        if os.path.exists(self._getPath('summary.txt')):
            fh=open(self._getPath('summary.txt'))
            for line in fh.readlines():
                summary.append(line.strip())
            fh.close()
        return summary

    def _citations(self):
        return ['Beckers2019']
