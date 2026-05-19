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
    """
    Applies statistical segmentation to three-dimensional cryo-EM density maps using
    False Discovery Rate and Family-Wise Error Rate control approaches. The protocol
    is designed to distinguish meaningful signal from background noise in a rigorous
    and reproducible manner, helping users identify density regions that are
    statistically supported while reducing the risk of interpreting random noise as
    structural information.

    AI Generated:

    FDR Segmentation (ProtFDR) - User Manual
        Overview

        The FDR Segmentation protocol applies statistical hypothesis testing to cryo-EM
        density maps in order to identify regions that are likely to represent true
        molecular signal. Rather than relying exclusively on visual inspection or
        manually selected thresholds, the protocol estimates significance levels from
        the noise characteristics of the map and generates a confidence-based
        segmentation.

        This approach is particularly valuable when working with weak densities,
        heterogeneous specimens, flexible regions, or maps where the distinction
        between signal and noise is not obvious. By controlling the expected rate of
        false positive detections, the protocol provides a statistically principled
        framework for map interpretation.

        Inputs and General Workflow

        The protocol requires a reconstructed cryo-EM volume as its primary input.
        Depending on the available data, users may analyze a single reconstructed map
        or make use of associated half maps when independent noise estimation is
        desired. The workflow evaluates the statistical properties of the density and
        produces a confidence map that highlights regions supported by the selected
        significance criterion.

        For biological interpretation, this process can help determine which features
        are sufficiently reliable for model building, structural analysis, or
        visualization. It is especially useful when assessing low-occupancy features,
        flexible domains, ligand densities, or peripheral regions that may be close
        to the noise level.

        Noise Estimation

        Accurate noise characterization is a critical component of the analysis. The
        protocol estimates background statistics from regions expected to contain
        solvent rather than molecular density. The user specifies a sampling region
        intended to represent noise, and the quality of the statistical assessment
        depends strongly on how representative this region is.

        From a biological perspective, the selected noise region should avoid any
        genuine molecular features and should reflect the background characteristics
        of the reconstruction. If the noise estimate is inaccurate, the resulting
        confidence values may become either overly permissive or excessively
        conservative.

        Use of Half Maps

        When independent half maps are available, they can be incorporated into the
        statistical estimation process. This option is often advantageous because it
        allows noise characteristics to be evaluated using independently reconstructed
        data, reducing bias and improving reliability.

        Half-map based analysis is particularly valuable in modern cryo-EM workflows,
        where independent reconstructions are routinely generated for validation
        purposes. In many cases, this approach provides greater confidence in the
        resulting segmentation and significance estimates.

        Choice of Statistical Method

        Several multiple-testing correction strategies are available. False Discovery
        Rate methods aim to control the expected proportion of false positive
        detections among accepted voxels, making them suitable for exploratory
        structural analyses where sensitivity is important.

        Family-Wise Error Rate methods apply stricter statistical control and are more
        conservative. These approaches reduce the probability of any false positive
        detection but may exclude weaker genuine signal. They are often preferred when
        particularly high confidence is required, such as in validation studies or
        when evaluating subtle structural features.

        The choice between these methods depends on the biological objective. Users
        seeking broader detection of potential features may favor FDR-based methods,
        whereas users prioritizing maximal statistical confidence may prefer FWER
        approaches.

        Outputs and Their Interpretation

        The principal result is a confidence map in which density values reflect the
        statistical support for the presence of signal. This output can be visualized
        alongside the original reconstruction and used to guide threshold selection,
        model building, or biological interpretation.

        The protocol also generates diagnostic information that helps evaluate the
        quality of the noise estimation and the resulting statistical analysis.
        Reviewing these diagnostics is strongly recommended before drawing biological
        conclusions.

        Practical Recommendations

        In routine cryo-EM analysis, users should carefully verify that the selected
        noise region is representative of background solvent and free of molecular
        density. Diagnostic outputs should always be inspected to confirm that the
        statistical assumptions are reasonable.

        For exploratory studies involving weak features or flexible regions, FDR-based
        methods often provide a useful balance between sensitivity and confidence. For
        highly conservative analyses, particularly when validating novel structural
        observations, FWER-based methods may be more appropriate.

        When half maps are available, their use is generally recommended because they
        provide a stronger foundation for statistical assessment and can improve the
        robustness of the resulting confidence estimates.

        Final Perspective

        Statistical segmentation provides an objective framework for interpreting
        cryo-EM density maps. By quantifying confidence rather than relying solely on
        visual thresholding, the protocol helps users distinguish reliable structural
        features from background fluctuations. Careful noise estimation, appropriate
        selection of the statistical correction method, and thoughtful biological
        interpretation are essential for obtaining meaningful and trustworthy results.
    """
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
