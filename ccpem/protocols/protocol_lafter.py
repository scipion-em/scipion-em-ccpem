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
    """
    Applies LAFTER as a post-processing filter to cryo-EM density maps. The protocol is designed to improve the interpretability of reconstructed volumes by reducing noise and enhancing meaningful structural signal while preserving information supported by independent half-map validation.

    AI Generated:

    LAFTER Post-Processing (ProtLafter) - User Manual
        Overview

        The LAFTER protocol performs post-processing of cryo-EM reconstructions using information derived from
        two independently reconstructed half maps. Its primary goal is to suppress noise while retaining
        structural features that are consistently supported by the experimental data. This approach can improve
        map quality and facilitate downstream interpretation, visualization, model building, and biological
        analysis.

        In cryo-EM workflows, post-processing is often one of the final stages before structural interpretation.
        At this point, the reconstruction may still contain residual noise, variable local signal quality,
        or artifacts resulting from overfitting. LAFTER provides an additional filtering strategy that can
        improve the visibility of biologically relevant features while maintaining consistency with the
        information present in the independent half-map reconstructions.

        Inputs and General Workflow

        The protocol requires a reconstructed volume associated with exactly two half maps. These half maps
        provide the independent information necessary to estimate signal reliability and distinguish genuine
        structural features from noise. Because the method relies on this comparison, volumes lacking valid
        half maps are not suitable for processing.

        During execution, the protocol analyzes the half-map information and generates filtered outputs that
        emphasize reproducible structural content. The resulting maps can then be used for visualization,
        interpretation, refinement, or model construction.

        Defining the Region of Interest

        A key aspect of successful post-processing is defining the molecular region that should guide the
        filtering procedure. The protocol supports either an explicit mask or a radius-based definition of
        the particle region.

        When a mask is provided, filtering is focused on the user-defined molecular volume. This option is
        generally preferred when a reliable mask already exists, particularly for irregularly shaped particles,
        flexible assemblies, membrane proteins, or complexes with substantial solvent regions.

        Alternatively, a particle radius can be used to approximate the molecular boundary. This approach is
        convenient when a dedicated mask is unavailable or when rapid exploratory processing is desired.
        Care should be taken to select a radius that encompasses the biologically relevant structure without
        including excessive solvent volume.

        Noise Suppression and Signal Preservation

        The principal objective of the protocol is to separate meaningful signal from noise. By leveraging
        agreement between independent half maps, the method attempts to preserve reproducible structural
        information while reducing features that are unlikely to correspond to genuine molecular density.

        This can be particularly useful when working with maps that contain weak features, heterogeneous
        regions, or residual reconstruction artifacts. Improved noise suppression often leads to cleaner
        density maps that are easier to interpret biologically and visually.

        Sharpening Options

        The protocol optionally supports sharpening during post-processing. Sharpening can enhance the
        visibility of high-resolution features such as secondary structure elements, side-chain density,
        and fine molecular details.

        Although sharpening may improve interpretability, excessive enhancement can also amplify noise or
        create misleading visual features. Biological users should therefore evaluate sharpened maps together
        with other validation metrics and compare them against the original reconstruction when necessary.

        FSC-Based Filtering

        Filtering behavior is guided by a Fourier Shell Correlation threshold that determines the estimated
        boundary between reliable signal and noise. The selected cutoff influences how aggressively the map
        is filtered and therefore affects the balance between feature preservation and noise reduction.

        In most situations, standard FSC thresholds provide appropriate results. Advanced users may adjust
        this parameter when working with unusual datasets or when specific validation strategies require
        alternative criteria.

        Overfitting Considerations

        The protocol includes an optional strategy aimed at reducing the impact of overfitting. This feature
        may be beneficial when reconstructions contain suspicious high-resolution detail that is not strongly
        supported by independent validation.

        From a biological perspective, reducing overfitting helps ensure that observed density features are
        more likely to represent genuine structural information rather than processing artifacts. Users should
        still evaluate results in conjunction with established validation procedures and visual inspection.

        Sampling and Resolution Considerations

        The protocol can optionally adjust the final sampling characteristics of the processed map. This
        flexibility allows users to generate outputs that either preserve intermediate processing resolution
        or match the sampling characteristics of the original reconstruction.

        The choice depends on downstream requirements. For visualization and comparison with existing data,
        maintaining consistency with the original sampling is often advantageous. For specialized analyses,
        alternative sampling behavior may sometimes be preferred.

        Outputs and Their Interpretation

        The protocol produces a filtered volume representing the primary post-processed reconstruction. This
        map is intended to provide improved signal-to-noise characteristics and enhanced interpretability for
        structural analysis.

        In addition, a noise-suppressed volume is generated. This output can be useful for evaluating the
        effects of filtering and for assessing how noise reduction influences the appearance of structural
        features.

        Practical Recommendations

        For most cryo-EM projects, it is advisable to begin with a well-defined mask whenever one is
        available. Accurate masking typically improves the reliability of post-processing and reduces the
        influence of solvent noise.

        Sharpening should be applied conservatively and evaluated visually. If high-resolution features
        become unusually prominent or inconsistent with known structural characteristics, additional
        validation is recommended.

        When uncertainty exists regarding potential overfitting, generating results with and without
        overfitting reduction can help determine whether observed features remain stable under different
        processing conditions.

        Final Perspective

        Effective post-processing is an important step between reconstruction and biological interpretation.
        By combining half-map validation, noise suppression, optional sharpening, and overfitting control,
        the LAFTER protocol provides a framework for generating cleaner and more interpretable cryo-EM maps.
        Careful selection of masking strategy and filtering parameters can substantially improve confidence
        in the structural conclusions drawn from the final reconstruction.
    """
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
