from oncomind.annotator.clients.myvariant import AnnotatorMyVariantClient
from oncomind.annotator.clients.fda import AnnotatorFDAClient
from oncomind.annotator.clients.cgi import AnnotatorCGIClient
from oncomind.annotator.clients.oncotree import AnnotatorOncoTreeClient
from oncomind.annotator.clients.vicc import AnnotatorVICCClient
from oncomind.annotator.clients.civic import AnnotatorCIViCClient
from oncomind.annotator.clients.cbioportal import AnnotatorCBioPortalClient
from oncomind.annotator.clients.depmap import AnnotatorDepMapClient
from oncomind.annotator.clients.hotspots import AnnotatorHotspotsClient
from oncomind.annotator.clients.clinicaltrials import AnnotatorClinicalTrialsClient
from oncomind.annotator.clients.pubmed import AnnotatorPubMedClient
from oncomind.annotator.clients.semantic_scholar import AnnotatorSemanticScholarClient

__all__ = [
    "AnnotatorMyVariantClient",
    "AnnotatorFDAClient",
    "AnnotatorCGIClient",
    "AnnotatorOncoTreeClient",
    "AnnotatorVICCClient",
    "AnnotatorCIViCClient",
    "AnnotatorCBioPortalClient",
    "AnnotatorDepMapClient",
    "AnnotatorHotspotsClient",
    "AnnotatorClinicalTrialsClient",
    "AnnotatorPubMedClient",
    "AnnotatorSemanticScholarClient",
]
