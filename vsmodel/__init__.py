"""vsmodel -- vowel-space model of Berthommier (JASA-EL / arXiv:2111.00868).

Modules
-------
tlm        transmission-line model of the vocal tract (P. Badin, ICP)
models     three-phase mixing, coordination function, DRM / Fant models
synthesis  LPC vowel resynthesis from TLM spectra (L. Girin, ICP)
"""

from . import tlm, models, synthesis

__all__ = ['tlm', 'models', 'synthesis']
__version__ = '1.0'
