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
from pyworkflow.protocol.params import (EnumParam, PointerParam, BooleanParam, IntParam, FloatParam)
from pwem.emlib.image import ImageHandler
from pwem.objects import Volume

class ProtLafter(EMProtocol):
    """Applies Lafter as a post-processing filter"""
    _label = 'lafter'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVol', PointerParam, pointerClass="Volume",
                       label='Volume to filter:', allowsNull=False,
                       help="It must have two half maps")
        form.addParam('maskType', EnumParam, choices=['Mask', 'Radius'], default=1,
                       label='Mask type:')
        form.addParam('volumeMask', PointerParam, pointerClass="VolumeMask",
                       label='Mask:', condition="maskType==0")
        form.addParam('volumeRadius', IntParam, default=-1,
                       label='Volume radius (px):', condition="maskType==1")
        form.addParam('sharp', BooleanParam, default=False,
                       label='Sharp', help='Do sharpening')
        form.addParam('fscCutoff', FloatParam, default=0.143, expertLevel=LEVEL_ADVANCED,
                       label='FSC Cutoff', help='Cutoff of the FSC')
        form.addParam('downsample', BooleanParam, default=True,
                       label='Downsample',
                       help='Do downsampling at the end of the process to have the same sampling rate as the input')
        form.addParam('overfitting', BooleanParam, default=False,
                       label='Overfitting', help='Try to reduce overfitting in the map')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('lafterStep')
        self._insertFunctionStep('createOutput')

    def lafterStep(self):
        halfmaps = self.inputVol.get().getHalfMaps().split(',')
        fn1=halfmaps[0]
        fn2=halfmaps[1]
        fn1mrc = self._getTmpPath('vol1.mrc')
        fn2mrc = self._getTmpPath('vol2.mrc')
        ih = ImageHandler()
        ih.convert(fn1,fn1mrc)
        ih.convert(fn2,fn2mrc)

        args='--v1 tmp/vol1.mrc --v2 tmp/vol2.mrc'
        if self.maskType.get()==0:
            fnMmrc = self._getTmpPath('mask.mrc')
            ih.convert(self.volumeMask.get().getFileName(),fnMmrc)
            args+=" --mask tmp/mask.mrc"
        else:
            if self.volumeRadius.get()>0:
                args+=" --particle_diameter %d"%(2*self.volumeRadius.get())
            else:
                args+=" --particle_diameter %d"%(self.inputVol.get().getDim()[0])
        if self.sharp.get():
            args+=" --sharp"
        args +=" --fsc %f"%self.fscCutoff.get()
        if self.downsample.get():
            args+=" --downsample"
        if self.overfitting.get():
            args+=" --overfitting"
        self.runJob('lafter',args,cwd=self._getPath())

    def createOutput(self):
        Ts = self.inputVol.get().getSamplingRate()
        if not self.downsample.get():
            Ts=Ts/2
        volumeLafter=Volume()
        volumeLafter.setFileName(self._getPath("LAFTER_filtered.mrc"))
        volumeLafter.setSamplingRate(Ts)
        self._defineOutputs(outputVolumeLafter=volumeLafter)
        self._defineSourceRelation(self.inputVol.get(),volumeLafter)

        volumeSuppressed=Volume()
        volumeSuppressed.setFileName(self._getPath("noise_suppressed.mrc"))
        volumeSuppressed.setSamplingRate(Ts)
        self._defineOutputs(outputVolumeSuppresed=volumeSuppressed)
        self._defineSourceRelation(self.inputVol.get(),volumeSuppressed)

    # --------------------------- UTILS functions ------------------
    def _validate(self):
        errors = []
        if not self.inputVol.hasValue():
            errors.append("You must provide an input volume")
            return errors
        halfmaps = self.inputVol.get().getHalfMaps().split(',')
        if len(halfmaps)!=2:
            errors.append("The input volume must have two associated half maps.")
        if self.maskType.get()==0 and not self.volumeMask.hasValue():
            errors.append("You must provide a mask")
        return errors

    def _summary(self):
        summary = []
        return summary

    def _citations(self):
        return ['Ramlaul2019']
