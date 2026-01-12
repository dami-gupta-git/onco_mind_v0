"""FDA Drug Lookup Module.

Provides utilities to look up FDA drug approval information from Drugs@FDA.

Primary source: https://www.accessdata.fda.gov/scripts/cder/daf/
Fallback: FDA press releases at https://www.fda.gov/news-events/press-announcements

Usage:
    from oncomind.api.fda_drugs import FDADrugLookup

    lookup = FDADrugLookup()
    info = lookup.get_drug_info("capivasertib")
    # Returns: {"drugs_at_fda_url": "...", "label_url": "...", "press_release_url": "..."}
"""

import logging
from typing import Any
from urllib.parse import quote_plus

logger = logging.getLogger(__name__)

# Curated mapping of oncology drugs to their FDA application numbers (NDA/BLA)
# Format: drug_name (lowercase) -> {"appl_no": str, "brand": str, "generic": str, "label_url": str | None}
# Application numbers allow direct linking to Drugs@FDA product pages
ONCOLOGY_DRUG_MAPPINGS: dict[str, dict[str, Any]] = {
    # AKT inhibitors
    "capivasertib": {
        "appl_no": "218197",
        "brand": "Truqap",
        "generic": "capivasertib",
        "approval_date": "2023-11-16",
    },
    "truqap": {
        "appl_no": "218197",
        "brand": "Truqap",
        "generic": "capivasertib",
        "approval_date": "2023-11-16",
    },
    "azd5363": {  # Development code for capivasertib
        "appl_no": "218197",
        "brand": "Truqap",
        "generic": "capivasertib",
        "approval_date": "2023-11-16",
    },

    # BRAF inhibitors
    "vemurafenib": {
        "appl_no": "202429",
        "brand": "Zelboraf",
        "generic": "vemurafenib",
        "approval_date": "2011-08-17",
    },
    "zelboraf": {
        "appl_no": "202429",
        "brand": "Zelboraf",
        "generic": "vemurafenib",
        "approval_date": "2011-08-17",
    },
    "dabrafenib": {
        "appl_no": "202806",
        "brand": "Tafinlar",
        "generic": "dabrafenib",
        "approval_date": "2013-05-29",
    },
    "tafinlar": {
        "appl_no": "202806",
        "brand": "Tafinlar",
        "generic": "dabrafenib",
        "approval_date": "2013-05-29",
    },
    "encorafenib": {
        "appl_no": "210496",
        "brand": "Braftovi",
        "generic": "encorafenib",
        "approval_date": "2018-06-27",
    },
    "braftovi": {
        "appl_no": "210496",
        "brand": "Braftovi",
        "generic": "encorafenib",
        "approval_date": "2018-06-27",
    },

    # MEK inhibitors
    "trametinib": {
        "appl_no": "204114",
        "brand": "Mekinist",
        "generic": "trametinib",
        "approval_date": "2013-05-29",
    },
    "mekinist": {
        "appl_no": "204114",
        "brand": "Mekinist",
        "generic": "trametinib",
        "approval_date": "2013-05-29",
    },
    "cobimetinib": {
        "appl_no": "206192",
        "brand": "Cotellic",
        "generic": "cobimetinib",
        "approval_date": "2015-11-10",
    },
    "cotellic": {
        "appl_no": "206192",
        "brand": "Cotellic",
        "generic": "cobimetinib",
        "approval_date": "2015-11-10",
    },
    "binimetinib": {
        "appl_no": "210498",
        "brand": "Mektovi",
        "generic": "binimetinib",
        "approval_date": "2018-06-27",
    },
    "mektovi": {
        "appl_no": "210498",
        "brand": "Mektovi",
        "generic": "binimetinib",
        "approval_date": "2018-06-27",
    },

    # EGFR inhibitors
    "erlotinib": {
        "appl_no": "021743",
        "brand": "Tarceva",
        "generic": "erlotinib",
        "approval_date": "2004-11-18",
    },
    "tarceva": {
        "appl_no": "021743",
        "brand": "Tarceva",
        "generic": "erlotinib",
        "approval_date": "2004-11-18",
    },
    "gefitinib": {
        "appl_no": "206995",
        "brand": "Iressa",
        "generic": "gefitinib",
        "approval_date": "2015-07-13",
    },
    "iressa": {
        "appl_no": "206995",
        "brand": "Iressa",
        "generic": "gefitinib",
        "approval_date": "2015-07-13",
    },
    "afatinib": {
        "appl_no": "201292",
        "brand": "Gilotrif",
        "generic": "afatinib",
        "approval_date": "2013-07-12",
    },
    "gilotrif": {
        "appl_no": "201292",
        "brand": "Gilotrif",
        "generic": "afatinib",
        "approval_date": "2013-07-12",
    },
    "osimertinib": {
        "appl_no": "208065",
        "brand": "Tagrisso",
        "generic": "osimertinib",
        "approval_date": "2015-11-13",
    },
    "tagrisso": {
        "appl_no": "208065",
        "brand": "Tagrisso",
        "generic": "osimertinib",
        "approval_date": "2015-11-13",
    },
    "azd9291": {  # Development code for osimertinib
        "appl_no": "208065",
        "brand": "Tagrisso",
        "generic": "osimertinib",
        "approval_date": "2015-11-13",
    },

    # ALK inhibitors
    "crizotinib": {
        "appl_no": "202570",
        "brand": "Xalkori",
        "generic": "crizotinib",
        "approval_date": "2011-08-26",
    },
    "xalkori": {
        "appl_no": "202570",
        "brand": "Xalkori",
        "generic": "crizotinib",
        "approval_date": "2011-08-26",
    },
    "ceritinib": {
        "appl_no": "205755",
        "brand": "Zykadia",
        "generic": "ceritinib",
        "approval_date": "2014-04-29",
    },
    "zykadia": {
        "appl_no": "205755",
        "brand": "Zykadia",
        "generic": "ceritinib",
        "approval_date": "2014-04-29",
    },
    "alectinib": {
        "appl_no": "208434",
        "brand": "Alecensa",
        "generic": "alectinib",
        "approval_date": "2015-12-11",
    },
    "alecensa": {
        "appl_no": "208434",
        "brand": "Alecensa",
        "generic": "alectinib",
        "approval_date": "2015-12-11",
    },
    "brigatinib": {
        "appl_no": "208772",
        "brand": "Alunbrig",
        "generic": "brigatinib",
        "approval_date": "2017-04-28",
    },
    "alunbrig": {
        "appl_no": "208772",
        "brand": "Alunbrig",
        "generic": "brigatinib",
        "approval_date": "2017-04-28",
    },
    "lorlatinib": {
        "appl_no": "210868",
        "brand": "Lorbrena",
        "generic": "lorlatinib",
        "approval_date": "2018-11-02",
    },
    "lorbrena": {
        "appl_no": "210868",
        "brand": "Lorbrena",
        "generic": "lorlatinib",
        "approval_date": "2018-11-02",
    },

    # KRAS G12C inhibitors
    "sotorasib": {
        "appl_no": "214665",
        "brand": "Lumakras",
        "generic": "sotorasib",
        "approval_date": "2021-05-28",
    },
    "lumakras": {
        "appl_no": "214665",
        "brand": "Lumakras",
        "generic": "sotorasib",
        "approval_date": "2021-05-28",
    },
    "amg 510": {  # Development code for sotorasib
        "appl_no": "214665",
        "brand": "Lumakras",
        "generic": "sotorasib",
        "approval_date": "2021-05-28",
    },
    "amg510": {  # Development code for sotorasib (no space)
        "appl_no": "214665",
        "brand": "Lumakras",
        "generic": "sotorasib",
        "approval_date": "2021-05-28",
    },
    "adagrasib": {
        "appl_no": "216340",
        "brand": "Krazati",
        "generic": "adagrasib",
        "approval_date": "2022-12-12",
    },
    "krazati": {
        "appl_no": "216340",
        "brand": "Krazati",
        "generic": "adagrasib",
        "approval_date": "2022-12-12",
    },
    "mrtx849": {  # Development code for adagrasib
        "appl_no": "216340",
        "brand": "Krazati",
        "generic": "adagrasib",
        "approval_date": "2022-12-12",
    },

    # HER2 inhibitors
    "trastuzumab": {
        "appl_no": "103792",
        "brand": "Herceptin",
        "generic": "trastuzumab",
        "approval_date": "1998-09-25",
    },
    "herceptin": {
        "appl_no": "103792",
        "brand": "Herceptin",
        "generic": "trastuzumab",
        "approval_date": "1998-09-25",
    },
    "pertuzumab": {
        "appl_no": "125409",
        "brand": "Perjeta",
        "generic": "pertuzumab",
        "approval_date": "2012-06-08",
    },
    "perjeta": {
        "appl_no": "125409",
        "brand": "Perjeta",
        "generic": "pertuzumab",
        "approval_date": "2012-06-08",
    },
    "lapatinib": {
        "appl_no": "022059",
        "brand": "Tykerb",
        "generic": "lapatinib",
        "approval_date": "2007-03-13",
    },
    "tykerb": {
        "appl_no": "022059",
        "brand": "Tykerb",
        "generic": "lapatinib",
        "approval_date": "2007-03-13",
    },
    "neratinib": {
        "appl_no": "208051",
        "brand": "Nerlynx",
        "generic": "neratinib",
        "approval_date": "2017-07-17",
    },
    "nerlynx": {
        "appl_no": "208051",
        "brand": "Nerlynx",
        "generic": "neratinib",
        "approval_date": "2017-07-17",
    },
    "tucatinib": {
        "appl_no": "213411",
        "brand": "Tukysa",
        "generic": "tucatinib",
        "approval_date": "2020-04-17",
    },
    "tukysa": {
        "appl_no": "213411",
        "brand": "Tukysa",
        "generic": "tucatinib",
        "approval_date": "2020-04-17",
    },
    "trastuzumab deruxtecan": {
        "appl_no": "761139",
        "brand": "Enhertu",
        "generic": "trastuzumab deruxtecan",
        "approval_date": "2019-12-20",
    },
    "enhertu": {
        "appl_no": "761139",
        "brand": "Enhertu",
        "generic": "trastuzumab deruxtecan",
        "approval_date": "2019-12-20",
    },

    # PARP inhibitors
    "olaparib": {
        "appl_no": "206162",
        "brand": "Lynparza",
        "generic": "olaparib",
        "approval_date": "2014-12-19",
    },
    "lynparza": {
        "appl_no": "206162",
        "brand": "Lynparza",
        "generic": "olaparib",
        "approval_date": "2014-12-19",
    },
    "rucaparib": {
        "appl_no": "209115",
        "brand": "Rubraca",
        "generic": "rucaparib",
        "approval_date": "2016-12-19",
    },
    "rubraca": {
        "appl_no": "209115",
        "brand": "Rubraca",
        "generic": "rucaparib",
        "approval_date": "2016-12-19",
    },
    "niraparib": {
        "appl_no": "208447",
        "brand": "Zejula",
        "generic": "niraparib",
        "approval_date": "2017-03-27",
    },
    "zejula": {
        "appl_no": "208447",
        "brand": "Zejula",
        "generic": "niraparib",
        "approval_date": "2017-03-27",
    },
    "talazoparib": {
        "appl_no": "211651",
        "brand": "Talzenna",
        "generic": "talazoparib",
        "approval_date": "2018-10-16",
    },
    "talzenna": {
        "appl_no": "211651",
        "brand": "Talzenna",
        "generic": "talazoparib",
        "approval_date": "2018-10-16",
    },

    # Checkpoint inhibitors
    "pembrolizumab": {
        "appl_no": "125514",
        "brand": "Keytruda",
        "generic": "pembrolizumab",
        "approval_date": "2014-09-04",
    },
    "keytruda": {
        "appl_no": "125514",
        "brand": "Keytruda",
        "generic": "pembrolizumab",
        "approval_date": "2014-09-04",
    },
    "nivolumab": {
        "appl_no": "125554",
        "brand": "Opdivo",
        "generic": "nivolumab",
        "approval_date": "2014-12-22",
    },
    "opdivo": {
        "appl_no": "125554",
        "brand": "Opdivo",
        "generic": "nivolumab",
        "approval_date": "2014-12-22",
    },
    "atezolizumab": {
        "appl_no": "761034",
        "brand": "Tecentriq",
        "generic": "atezolizumab",
        "approval_date": "2016-05-18",
    },
    "tecentriq": {
        "appl_no": "761034",
        "brand": "Tecentriq",
        "generic": "atezolizumab",
        "approval_date": "2016-05-18",
    },
    "durvalumab": {
        "appl_no": "761069",
        "brand": "Imfinzi",
        "generic": "durvalumab",
        "approval_date": "2017-05-01",
    },
    "imfinzi": {
        "appl_no": "761069",
        "brand": "Imfinzi",
        "generic": "durvalumab",
        "approval_date": "2017-05-01",
    },
    "ipilimumab": {
        "appl_no": "125377",
        "brand": "Yervoy",
        "generic": "ipilimumab",
        "approval_date": "2011-03-25",
    },
    "yervoy": {
        "appl_no": "125377",
        "brand": "Yervoy",
        "generic": "ipilimumab",
        "approval_date": "2011-03-25",
    },

    # BCR-ABL inhibitors
    "imatinib": {
        "appl_no": "021588",
        "brand": "Gleevec",
        "generic": "imatinib",
        "approval_date": "2001-05-10",
    },
    "gleevec": {
        "appl_no": "021588",
        "brand": "Gleevec",
        "generic": "imatinib",
        "approval_date": "2001-05-10",
    },
    "dasatinib": {
        "appl_no": "021986",
        "brand": "Sprycel",
        "generic": "dasatinib",
        "approval_date": "2006-06-28",
    },
    "sprycel": {
        "appl_no": "021986",
        "brand": "Sprycel",
        "generic": "dasatinib",
        "approval_date": "2006-06-28",
    },
    "nilotinib": {
        "appl_no": "022068",
        "brand": "Tasigna",
        "generic": "nilotinib",
        "approval_date": "2007-10-29",
    },
    "tasigna": {
        "appl_no": "022068",
        "brand": "Tasigna",
        "generic": "nilotinib",
        "approval_date": "2007-10-29",
    },
    "bosutinib": {
        "appl_no": "203341",
        "brand": "Bosulif",
        "generic": "bosutinib",
        "approval_date": "2012-09-04",
    },
    "bosulif": {
        "appl_no": "203341",
        "brand": "Bosulif",
        "generic": "bosutinib",
        "approval_date": "2012-09-04",
    },
    "ponatinib": {
        "appl_no": "203469",
        "brand": "Iclusig",
        "generic": "ponatinib",
        "approval_date": "2012-12-14",
    },
    "iclusig": {
        "appl_no": "203469",
        "brand": "Iclusig",
        "generic": "ponatinib",
        "approval_date": "2012-12-14",
    },

    # CDK4/6 inhibitors
    "palbociclib": {
        "appl_no": "207103",
        "brand": "Ibrance",
        "generic": "palbociclib",
        "approval_date": "2015-02-03",
    },
    "ibrance": {
        "appl_no": "207103",
        "brand": "Ibrance",
        "generic": "palbociclib",
        "approval_date": "2015-02-03",
    },
    "ribociclib": {
        "appl_no": "209092",
        "brand": "Kisqali",
        "generic": "ribociclib",
        "approval_date": "2017-03-13",
    },
    "kisqali": {
        "appl_no": "209092",
        "brand": "Kisqali",
        "generic": "ribociclib",
        "approval_date": "2017-03-13",
    },
    "abemaciclib": {
        "appl_no": "208716",
        "brand": "Verzenio",
        "generic": "abemaciclib",
        "approval_date": "2017-09-28",
    },
    "verzenio": {
        "appl_no": "208716",
        "brand": "Verzenio",
        "generic": "abemaciclib",
        "approval_date": "2017-09-28",
    },

    # PI3K inhibitors
    "alpelisib": {
        "appl_no": "212526",
        "brand": "Piqray",
        "generic": "alpelisib",
        "approval_date": "2019-05-24",
    },
    "piqray": {
        "appl_no": "212526",
        "brand": "Piqray",
        "generic": "alpelisib",
        "approval_date": "2019-05-24",
    },

    # RET inhibitors
    "selpercatinib": {
        "appl_no": "213246",
        "brand": "Retevmo",
        "generic": "selpercatinib",
        "approval_date": "2020-05-08",
    },
    "retevmo": {
        "appl_no": "213246",
        "brand": "Retevmo",
        "generic": "selpercatinib",
        "approval_date": "2020-05-08",
    },
    "pralsetinib": {
        "appl_no": "214701",
        "brand": "Gavreto",
        "generic": "pralsetinib",
        "approval_date": "2020-09-04",
    },
    "gavreto": {
        "appl_no": "214701",
        "brand": "Gavreto",
        "generic": "pralsetinib",
        "approval_date": "2020-09-04",
    },

    # MET inhibitors
    "capmatinib": {
        "appl_no": "213591",
        "brand": "Tabrecta",
        "generic": "capmatinib",
        "approval_date": "2020-05-06",
    },
    "tabrecta": {
        "appl_no": "213591",
        "brand": "Tabrecta",
        "generic": "capmatinib",
        "approval_date": "2020-05-06",
    },
    "tepotinib": {
        "appl_no": "214096",
        "brand": "Tepmetko",
        "generic": "tepotinib",
        "approval_date": "2021-02-03",
    },
    "tepmetko": {
        "appl_no": "214096",
        "brand": "Tepmetko",
        "generic": "tepotinib",
        "approval_date": "2021-02-03",
    },

    # EZH2 inhibitors
    "tazemetostat": {
        "appl_no": "213400",
        "brand": "Tazverik",
        "generic": "tazemetostat",
        "approval_date": "2020-01-23",
    },
    "tazverik": {
        "appl_no": "213400",
        "brand": "Tazverik",
        "generic": "tazemetostat",
        "approval_date": "2020-01-23",
    },

    # IDH inhibitors
    "ivosidenib": {
        "appl_no": "211192",
        "brand": "Tibsovo",
        "generic": "ivosidenib",
        "approval_date": "2018-07-20",
    },
    "tibsovo": {
        "appl_no": "211192",
        "brand": "Tibsovo",
        "generic": "ivosidenib",
        "approval_date": "2018-07-20",
    },
    "enasidenib": {
        "appl_no": "209606",
        "brand": "Idhifa",
        "generic": "enasidenib",
        "approval_date": "2017-08-01",
    },
    "idhifa": {
        "appl_no": "209606",
        "brand": "Idhifa",
        "generic": "enasidenib",
        "approval_date": "2017-08-01",
    },

    # FGFR inhibitors
    "erdafitinib": {
        "appl_no": "212018",
        "brand": "Balversa",
        "generic": "erdafitinib",
        "approval_date": "2019-04-12",
    },
    "balversa": {
        "appl_no": "212018",
        "brand": "Balversa",
        "generic": "erdafitinib",
        "approval_date": "2019-04-12",
    },
    "pemigatinib": {
        "appl_no": "213736",
        "brand": "Pemazyre",
        "generic": "pemigatinib",
        "approval_date": "2020-04-17",
    },
    "pemazyre": {
        "appl_no": "213736",
        "brand": "Pemazyre",
        "generic": "pemigatinib",
        "approval_date": "2020-04-17",
    },

    # NTRK inhibitors
    "larotrectinib": {
        "appl_no": "211710",
        "brand": "Vitrakvi",
        "generic": "larotrectinib",
        "approval_date": "2018-11-26",
    },
    "vitrakvi": {
        "appl_no": "211710",
        "brand": "Vitrakvi",
        "generic": "larotrectinib",
        "approval_date": "2018-11-26",
    },
    "entrectinib": {
        "appl_no": "212725",
        "brand": "Rozlytrek",
        "generic": "entrectinib",
        "approval_date": "2019-08-15",
    },
    "rozlytrek": {
        "appl_no": "212725",
        "brand": "Rozlytrek",
        "generic": "entrectinib",
        "approval_date": "2019-08-15",
    },

    # Fulvestrant (ER antagonist)
    "fulvestrant": {
        "appl_no": "021344",
        "brand": "Faslodex",
        "generic": "fulvestrant",
        "approval_date": "2002-04-25",
    },
    "faslodex": {
        "appl_no": "021344",
        "brand": "Faslodex",
        "generic": "fulvestrant",
        "approval_date": "2002-04-25",
    },
}


class FDADrugLookup:
    """Utility class for FDA drug information lookup."""

    # URL templates
    DRUGS_AT_FDA_PRODUCT_URL = "https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=overview.process&ApplNo={appl_no}"
    DRUGS_AT_FDA_SEARCH_URL = "https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&searchTerm={drug_name}"
    FDA_PRESS_SEARCH_URL = "https://www.fda.gov/news-events/press-announcements?search={drug_name}"
    FDA_SITE_SEARCH_URL = "https://search.usa.gov/search?affiliate=fda1&query={drug_name}"

    def __init__(self):
        self.mappings = ONCOLOGY_DRUG_MAPPINGS

    def normalize_drug_name(self, drug_name: str) -> str:
        """Normalize drug name for lookup."""
        # Lowercase and strip whitespace
        normalized = drug_name.lower().strip()
        # Remove common suffixes/prefixes
        for suffix in [" hydrochloride", " hcl", " mesylate", " tosylate", " maleate"]:
            if normalized.endswith(suffix):
                normalized = normalized[:-len(suffix)]
        return normalized

    def get_drug_info(self, drug_name: str) -> dict[str, Any]:
        """Get FDA drug information for a given drug name.

        Args:
            drug_name: Brand name, generic name, or code name

        Returns:
            Dictionary with:
                - found: bool - whether drug was found in curated mapping
                - brand: str | None - brand name
                - generic: str | None - generic name
                - appl_no: str | None - FDA application number
                - approval_date: str | None - approval date
                - drugs_at_fda_url: str - direct URL to Drugs@FDA page (if found) or search URL
                - drugs_at_fda_search_url: str - search URL for Drugs@FDA
                - fda_press_search_url: str - FDA press release search URL
                - fda_site_search_url: str - FDA site-wide search URL
        """
        normalized = self.normalize_drug_name(drug_name)

        # Check curated mapping
        if normalized in self.mappings:
            mapping = self.mappings[normalized]
            appl_no = mapping["appl_no"]
            return {
                "found": True,
                "brand": mapping.get("brand"),
                "generic": mapping.get("generic"),
                "appl_no": appl_no,
                "approval_date": mapping.get("approval_date"),
                "drugs_at_fda_url": self.DRUGS_AT_FDA_PRODUCT_URL.format(appl_no=appl_no),
                "drugs_at_fda_search_url": self.DRUGS_AT_FDA_SEARCH_URL.format(drug_name=quote_plus(drug_name)),
                "fda_press_search_url": self.FDA_PRESS_SEARCH_URL.format(drug_name=quote_plus(drug_name)),
                "fda_site_search_url": self.FDA_SITE_SEARCH_URL.format(drug_name=quote_plus(drug_name)),
            }

        # Not in curated mapping - return search URLs
        return {
            "found": False,
            "brand": None,
            "generic": None,
            "appl_no": None,
            "approval_date": None,
            "drugs_at_fda_url": None,  # Can't generate direct URL without appl_no
            "drugs_at_fda_search_url": self.DRUGS_AT_FDA_SEARCH_URL.format(drug_name=quote_plus(drug_name)),
            "fda_press_search_url": self.FDA_PRESS_SEARCH_URL.format(drug_name=quote_plus(drug_name)),
            "fda_site_search_url": self.FDA_SITE_SEARCH_URL.format(drug_name=quote_plus(drug_name)),
        }

    def get_drugs_at_fda_url(self, drug_name: str) -> str | None:
        """Get direct Drugs@FDA URL if drug is in curated mapping.

        Returns:
            Direct URL to product page, or None if not found
        """
        info = self.get_drug_info(drug_name)
        return info.get("drugs_at_fda_url")

    def get_best_fda_url(self, drug_name: str) -> str:
        """Get the best available FDA URL for a drug.

        Returns direct Drugs@FDA URL if available, otherwise search URL.
        """
        info = self.get_drug_info(drug_name)
        if info["drugs_at_fda_url"]:
            return info["drugs_at_fda_url"]
        return info["drugs_at_fda_search_url"]

    def is_known_drug(self, drug_name: str) -> bool:
        """Check if drug is in the curated mapping."""
        normalized = self.normalize_drug_name(drug_name)
        return normalized in self.mappings


# Module-level singleton for convenience
_lookup_instance: FDADrugLookup | None = None


def get_fda_drug_lookup() -> FDADrugLookup:
    """Get the singleton FDADrugLookup instance."""
    global _lookup_instance
    if _lookup_instance is None:
        _lookup_instance = FDADrugLookup()
    return _lookup_instance


def get_drugs_at_fda_url(drug_name: str) -> str | None:
    """Convenience function to get Drugs@FDA URL for a drug."""
    return get_fda_drug_lookup().get_drugs_at_fda_url(drug_name)


def get_best_fda_url(drug_name: str) -> str:
    """Convenience function to get best FDA URL for a drug."""
    return get_fda_drug_lookup().get_best_fda_url(drug_name)


# =============================================================================
# OpenFDA Label API Client
# =============================================================================

import re
import httpx
from pydantic import BaseModel, Field


class FDALabelInfo(BaseModel):
    """FDA drug label information from OpenFDA API.

    This provides authoritative, up-to-date FDA label data for cross-source validation.
    """

    # Drug identification
    drug_name: str
    brand_name: str | None = None
    generic_name: str | None = None

    # FDA identifiers
    application_number: str | None = Field(
        default=None,
        description="FDA application number (NDA/BLA/ANDA)"
    )
    set_id: str | None = Field(
        default=None,
        description="Unique identifier for this label set"
    )

    # URLs
    drugs_at_fda_url: str | None = Field(
        default=None,
        description="Direct URL to Drugs@FDA product page"
    )
    label_pdf_url: str | None = Field(
        default=None,
        description="URL to FDA label PDF (from DailyMed)"
    )
    pubmed_url: str | None = Field(
        default=None,
        description="PubMed search URL for the drug"
    )
    approval_docs: list[str] = Field(
        default_factory=list,
        description="URLs to approval letters, medical reviews, etc."
    )

    # Approval information
    approval_date: str | None = None
    approval_date_source: str = Field(
        default="openfda",
        description="Source of approval date: 'openfda', 'curated', or 'label_text'"
    )

    # Label content
    indications_and_usage: str | None = Field(
        default=None,
        description="Full indications and usage text from label"
    )
    biomarker_text_exact: str | None = Field(
        default=None,
        description="Extracted biomarker-specific indication text"
    )

    # Biomarker specificity - is approval for gene-level or specific variant?
    biomarker_specificity: str | None = Field(
        default=None,
        description="'variant' if approval is for specific variant (e.g., BRAF V600E), "
                    "'gene' if for any mutation in gene (e.g., BRCA-mutated), "
                    "'phenotype' if for phenotype (e.g., MSI-H, HRD)"
    )
    specific_variants: list[str] = Field(
        default_factory=list,
        description="List of specific variants mentioned in approval (e.g., ['V600E', 'V600K'])"
    )
    target_genes: list[str] = Field(
        default_factory=list,
        description="List of genes mentioned in approval"
    )

    # Extracted clinical data
    trials_cited: list[str] = Field(
        default_factory=list,
        description="NCT IDs of trials cited in the label"
    )

    # Metadata
    effective_date: str | None = Field(
        default=None,
        description="Effective date of this label version"
    )
    api_success: bool = Field(
        default=False,
        description="Whether OpenFDA API call succeeded"
    )
    error_message: str | None = None


class OpenFDAClient:
    """Client for OpenFDA Drug Label API.

    Documentation: https://open.fda.gov/apis/drug/label/
    """

    BASE_URL = "https://api.fda.gov/drug/label.json"
    DAILYMED_LABEL_URL = "https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid={set_id}"
    PUBMED_SEARCH_URL = "https://pubmed.ncbi.nlm.nih.gov/?term={drug_name}+cancer"

    # Patterns for extracting biomarker text
    BIOMARKER_PATTERNS = [
        # Gene + variant patterns
        r'((?:EGFR|BRAF|KRAS|NRAS|HER2|ERBB2|ALK|ROS1|RET|MET|NTRK|BRCA|PIK3CA|AKT1|IDH1|IDH2|FGFR|KIT|PDGFRA|FLT3|NPM1|TP53|ATM|PALB2|RAD51|CHEK2|CDH1|STK11|NF1|PTEN|RB1|SMAD4|ARID1A|MSH2|MLH1|MSH6|PMS2|POLE|POLD1)[^.]*(?:mutation|mutant|mutated|alteration|amplification|fusion|rearrangement|deletion|insertion|exon)[^.]*\.)',
        # MSI-H/dMMR patterns
        r'(microsatellite instability[^.]*\.)',
        r'(mismatch repair deficien[^.]*\.)',
        r'(MSI-H[^.]*\.)',
        r'(dMMR[^.]*\.)',
        # HRD patterns
        r'(homologous recombination[^.]*\.)',
        r'(HRD[^.]*\.)',
        # TMB patterns
        r'(tumor mutational burden[^.]*\.)',
        r'(TMB[^.]*\.)',
    ]

    # Pattern for NCT IDs
    NCT_PATTERN = r'NCT\d{8}'

    def __init__(self, timeout: float = 10.0):
        self.timeout = timeout
        self._drug_lookup = get_fda_drug_lookup()

    async def get_label_info(self, drug_name: str) -> FDALabelInfo:
        """Fetch FDA label information for a drug.

        Tries multiple search strategies:
        1. Generic name search
        2. Brand name search
        3. Substance name search

        Args:
            drug_name: Drug name (generic or brand)

        Returns:
            FDALabelInfo with available data
        """
        result = FDALabelInfo(drug_name=drug_name)

        # Try to get info from curated mapping first
        curated_info = self._drug_lookup.get_drug_info(drug_name)
        if curated_info["found"]:
            result.brand_name = curated_info.get("brand")
            result.generic_name = curated_info.get("generic")
            result.application_number = curated_info.get("appl_no")
            result.drugs_at_fda_url = curated_info.get("drugs_at_fda_url")
            if curated_info.get("approval_date"):
                result.approval_date = curated_info["approval_date"]
                result.approval_date_source = "curated"

        # Try OpenFDA API
        label_data = await self._fetch_label(drug_name)

        if label_data:
            result.api_success = True
            self._parse_label_data(label_data, result)
        else:
            # If first search failed and we have curated info, try with those names
            if curated_info["found"]:
                for name in [curated_info.get("generic"), curated_info.get("brand")]:
                    if name and name.lower() != drug_name.lower():
                        label_data = await self._fetch_label(name)
                        if label_data:
                            result.api_success = True
                            self._parse_label_data(label_data, result)
                            break

        if not result.api_success and not curated_info["found"]:
            result.error_message = f"Drug '{drug_name}' not found in OpenFDA or curated mapping"

        # Generate PubMed search URL using generic name if available
        search_term = result.generic_name or drug_name
        from urllib.parse import quote
        result.pubmed_url = self.PUBMED_SEARCH_URL.format(drug_name=quote(search_term))

        return result

    async def _fetch_label(self, drug_name: str) -> dict | None:
        """Fetch label data from OpenFDA API."""
        normalized = drug_name.lower().strip()

        # Try different search strategies
        search_queries = [
            f'openfda.generic_name:"{normalized}"',
            f'openfda.brand_name:"{normalized}"',
            f'openfda.substance_name:"{normalized}"',
        ]

        async with httpx.AsyncClient(timeout=self.timeout) as client:
            for query in search_queries:
                try:
                    response = await client.get(
                        self.BASE_URL,
                        params={"search": query, "limit": 1}
                    )

                    if response.status_code == 200:
                        data = response.json()
                        if data.get("results"):
                            return data["results"][0]

                except httpx.TimeoutException:
                    logger.warning(f"OpenFDA timeout for query: {query}")
                except httpx.HTTPError as e:
                    logger.warning(f"OpenFDA HTTP error: {e}")
                except Exception as e:
                    logger.warning(f"OpenFDA error: {e}")

        return None

    def _parse_label_data(self, data: dict, result: FDALabelInfo) -> None:
        """Parse OpenFDA label response into FDALabelInfo."""
        openfda = data.get("openfda", {})

        # Drug names
        if not result.brand_name and openfda.get("brand_name"):
            result.brand_name = openfda["brand_name"][0]
        if not result.generic_name and openfda.get("generic_name"):
            result.generic_name = openfda["generic_name"][0]

        # Application number
        if openfda.get("application_number"):
            appl_no = openfda["application_number"][0]
            result.application_number = appl_no
            # Build Drugs@FDA URL if we don't have one
            if not result.drugs_at_fda_url:
                # Extract just the number part (remove NDA/BLA prefix)
                num_match = re.search(r'\d+', appl_no)
                if num_match:
                    result.drugs_at_fda_url = self._drug_lookup.DRUGS_AT_FDA_PRODUCT_URL.format(
                        appl_no=num_match.group()
                    )

        # Set ID for DailyMed label URL
        if data.get("set_id"):
            result.set_id = data["set_id"]
            result.label_pdf_url = self.DAILYMED_LABEL_URL.format(set_id=data["set_id"])

        # Effective date
        if data.get("effective_time"):
            result.effective_date = data["effective_time"]

        # Indications and usage
        if data.get("indications_and_usage"):
            indications = data["indications_and_usage"]
            if isinstance(indications, list):
                result.indications_and_usage = " ".join(indications)
            else:
                result.indications_and_usage = str(indications)

            # Extract biomarker-specific text
            result.biomarker_text_exact = self._extract_biomarker_text(
                result.indications_and_usage
            )

            # Extract biomarker specificity (gene vs variant vs phenotype level)
            self._apply_biomarker_specificity(result.indications_and_usage, result)

        # Extract NCT IDs from clinical studies section
        clinical_studies = data.get("clinical_studies", [])
        if clinical_studies:
            studies_text = " ".join(clinical_studies) if isinstance(clinical_studies, list) else str(clinical_studies)
            result.trials_cited = self._extract_nct_ids(studies_text)

        # Also check indications for NCT IDs
        if result.indications_and_usage:
            additional_ncts = self._extract_nct_ids(result.indications_and_usage)
            for nct in additional_ncts:
                if nct not in result.trials_cited:
                    result.trials_cited.append(nct)

    def _extract_biomarker_text(self, text: str) -> str | None:
        """Extract biomarker-specific sentences from indication text."""
        if not text:
            return None

        biomarker_sentences = []

        for pattern in self.BIOMARKER_PATTERNS:
            matches = re.findall(pattern, text, re.IGNORECASE)
            for match in matches:
                cleaned = match.strip()
                if cleaned and cleaned not in biomarker_sentences:
                    biomarker_sentences.append(cleaned)

        if biomarker_sentences:
            return " ".join(biomarker_sentences)

        return None

    def _extract_nct_ids(self, text: str) -> list[str]:
        """Extract NCT trial IDs from text."""
        if not text:
            return []

        matches = re.findall(self.NCT_PATTERN, text)
        # Deduplicate while preserving order
        seen = set()
        unique = []
        for nct in matches:
            if nct not in seen:
                seen.add(nct)
                unique.append(nct)

        return unique

    def _apply_biomarker_specificity(self, text: str, result: FDALabelInfo) -> None:
        """Apply biomarker specificity parsing to FDALabelInfo result.

        Uses the consolidated parse_biomarker_specificity function and maps
        the result to FDALabelInfo fields.
        """
        if not text:
            return

        from oncomind.insight_builder.fda_processor import parse_biomarker_specificity

        text_upper = text.upper()

        # Known oncology genes - scan text for all mentioned genes
        oncology_genes = [
            'EGFR', 'BRAF', 'KRAS', 'NRAS', 'HER2', 'ERBB2', 'ALK', 'ROS1', 'RET',
            'MET', 'NTRK', 'NTRK1', 'NTRK2', 'NTRK3', 'BRCA1', 'BRCA2', 'BRCA',
            'PIK3CA', 'AKT1', 'IDH1', 'IDH2', 'FGFR', 'FGFR1', 'FGFR2', 'FGFR3',
            'KIT', 'PDGFRA', 'FLT3', 'NPM1', 'TP53', 'ATM', 'PALB2', 'RAD51',
            'CHEK2', 'CDH1', 'STK11', 'NF1', 'PTEN', 'RB1', 'SMAD4', 'ARID1A',
            'MSH2', 'MLH1', 'MSH6', 'PMS2', 'POLE', 'POLD1'
        ]

        # Find all genes mentioned in text
        genes_found = set()
        for gene in oncology_genes:
            gene_pattern = rf'\b{gene}\b'
            if re.search(gene_pattern, text_upper):
                genes_found.add(gene)

        # Try phenotype parsing first (no gene required)
        spec = parse_biomarker_specificity(text, gene=None)
        if spec and spec["level"] == "phenotype":
            result.biomarker_specificity = "phenotype"
            result.specific_variants = spec.get("phenotypes", [])
            result.target_genes = list(genes_found)
            return

        # Try each gene found to get variant-level specificity
        for gene in genes_found:
            spec = parse_biomarker_specificity(text, gene=gene)
            if spec:
                if spec["level"] == "variant":
                    result.biomarker_specificity = "variant"
                    result.specific_variants = [spec["specified_variant"]]
                    result.target_genes = spec.get("target_genes", [gene])
                    return
                elif spec["level"] == "variant_list":
                    result.biomarker_specificity = "variant"
                    result.specific_variants = spec.get("specified_variants", [])
                    result.target_genes = spec.get("target_genes", [gene])
                    return
                elif spec["level"] in ("codon", "exon"):
                    # Treat codon/exon as variant-level for this model
                    result.biomarker_specificity = "variant"
                    if spec["level"] == "codon":
                        result.specific_variants = [spec["codon"]]
                    else:
                        result.specific_variants = [f"exon {spec['exon']}"]
                    result.target_genes = spec.get("target_genes", [gene])
                    return
                elif spec["level"] == "gene":
                    result.biomarker_specificity = "gene"
                    result.target_genes = spec.get("target_genes", [gene])
                    # Don't return - continue looking for more specific matches
                    continue
                elif spec["level"] in ("contraindication", "wild_type_required"):
                    # These indicate drug is NOT for mutated gene
                    result.biomarker_specificity = spec["level"]
                    result.target_genes = spec.get("target_genes", [gene])
                    return

        # If we found any genes but no specific match, mark as gene-level
        if genes_found and not result.biomarker_specificity:
            result.biomarker_specificity = "gene"

        result.target_genes = list(genes_found)

    async def search_by_biomarker(
        self,
        gene: str,
        variant: str | None = None,
        limit: int = 100,
    ) -> list[FDALabelInfo]:
        """Search FDA drug labels by biomarker (gene/variant).

        Searches the indications_and_usage field for drugs that mention
        the specified gene. This catches FDA approvals that may not be
        in curated databases yet.

        Args:
            gene: Gene symbol (e.g., "IDH1", "BRAF", "EGFR")
            variant: Optional variant (e.g., "R132H", "V600E") to filter results
            limit: Maximum number of results to return

        Returns:
            List of FDALabelInfo objects for drugs mentioning the gene

        Example:
            >>> client = OpenFDAClient()
            >>> labels = await client.search_by_biomarker("IDH1")
            >>> for label in labels:
            ...     print(f"{label.brand_name}: {label.target_genes}")
        """
        gene_upper = gene.upper()
        results: list[FDALabelInfo] = []

        # Search indications_and_usage for the gene
        # Use exact word boundary matching where possible
        search_query = f'indications_and_usage:"{gene_upper}"'

        async with httpx.AsyncClient(timeout=self.timeout) as client:
            try:
                response = await client.get(
                    self.BASE_URL,
                    params={"search": search_query, "limit": limit}
                )

                if response.status_code != 200:
                    logger.warning(f"OpenFDA biomarker search failed: {response.status_code}")
                    return results

                data = response.json()
                api_results = data.get("results", [])

                for label_data in api_results:
                    # Parse into FDALabelInfo
                    openfda = label_data.get("openfda", {})

                    brand_name = openfda.get("brand_name", [None])[0] if openfda.get("brand_name") else None
                    generic_name = openfda.get("generic_name", [None])[0] if openfda.get("generic_name") else None
                    drug_name = generic_name or brand_name or "Unknown"

                    result = FDALabelInfo(drug_name=drug_name)
                    result.brand_name = brand_name
                    result.generic_name = generic_name
                    result.api_success = True

                    # Parse full label data
                    self._parse_label_data(label_data, result)

                    # Apply biomarker specificity analysis
                    if result.indications_and_usage:
                        self._apply_biomarker_specificity(result.indications_and_usage, result)

                    # Verify the gene is actually in target_genes
                    # (API search may have false positives)
                    if gene_upper not in [g.upper() for g in result.target_genes]:
                        # Check if gene appears in indications text
                        if result.indications_and_usage:
                            text_upper = result.indications_and_usage.upper()
                            if gene_upper not in text_upper:
                                continue
                            # Gene is in text but not parsed - add it
                            result.target_genes.append(gene_upper)

                    # If variant specified, filter by variant match
                    if variant:
                        variant_upper = variant.upper()
                        # Check if variant is in specific_variants or indications text
                        variant_matched = False

                        if result.specific_variants:
                            for v in result.specific_variants:
                                if variant_upper in v.upper():
                                    variant_matched = True
                                    break

                        # Also check indications text for variant mention
                        if not variant_matched and result.indications_and_usage:
                            if variant_upper in result.indications_and_usage.upper():
                                variant_matched = True

                        # For gene-level approvals (any mutation), include them
                        if not variant_matched:
                            if result.biomarker_specificity == "gene":
                                # Gene-level approval covers all variants
                                variant_matched = True
                            elif "MUTATION" in (result.indications_and_usage or "").upper():
                                # Generic mutation language
                                variant_matched = True

                        if not variant_matched:
                            continue

                    results.append(result)

            except httpx.TimeoutException:
                logger.warning(f"OpenFDA biomarker search timeout for gene: {gene}")
            except httpx.HTTPError as e:
                logger.warning(f"OpenFDA biomarker search HTTP error: {e}")
            except Exception as e:
                logger.error(f"OpenFDA biomarker search error: {e}")

        # Deduplicate by generic_name or brand_name
        seen_drugs: set[str] = set()
        unique_results: list[FDALabelInfo] = []
        for result in results:
            key = (result.generic_name or result.brand_name or result.drug_name).lower()
            if key not in seen_drugs:
                seen_drugs.add(key)
                unique_results.append(result)

        logger.info(f"OpenFDA biomarker search for {gene}: found {len(unique_results)} drugs")
        return unique_results


async def search_fda_by_biomarker(
    gene: str,
    variant: str | None = None,
    limit: int = 100,
) -> list[FDALabelInfo]:
    """Convenience function to search FDA labels by biomarker.

    Args:
        gene: Gene symbol (e.g., "IDH1", "BRAF")
        variant: Optional variant to filter results
        limit: Maximum results

    Returns:
        List of FDALabelInfo for drugs targeting this biomarker
    """
    client = get_openfda_client()
    return await client.search_by_biomarker(gene, variant, limit)


# Module-level singleton for OpenFDA client
_openfda_client: OpenFDAClient | None = None


def get_openfda_client() -> OpenFDAClient:
    """Get the singleton OpenFDAClient instance."""
    global _openfda_client
    if _openfda_client is None:
        _openfda_client = OpenFDAClient()
    return _openfda_client


async def get_fda_label_info(drug_name: str) -> FDALabelInfo:
    """Convenience function to get FDA label info for a drug."""
    return await get_openfda_client().get_label_info(drug_name)


# =============================================================================
# FDA Label Cache (for nightly batch population)
# =============================================================================

import json
from pathlib import Path
from datetime import datetime

from oncomind.config.constants import FDA_LABELS_FILE

# Default cache location
FDA_LABEL_CACHE_PATH = FDA_LABELS_FILE


class FDALabelCache:
    """File-based cache for FDA label information.

    Designed to be populated by a nightly batch job and read at runtime.
    This avoids API calls during user requests.
    """

    def __init__(self, cache_path: Path | None = None):
        self.cache_path = cache_path or FDA_LABEL_CACHE_PATH
        self._cache: dict[str, dict] = {}
        self._loaded = False
        self._last_updated: str | None = None

    def _ensure_loaded(self) -> None:
        """Load cache from disk if not already loaded."""
        if self._loaded:
            return

        if self.cache_path.exists():
            try:
                with open(self.cache_path, "r") as f:
                    data = json.load(f)
                    self._cache = data.get("drugs", {})
                    self._last_updated = data.get("last_updated")
                    logger.info(f"Loaded FDA label cache: {len(self._cache)} drugs, updated {self._last_updated}")
            except Exception as e:
                logger.warning(f"Failed to load FDA label cache: {e}")
                self._cache = {}

        self._loaded = True

    def get(self, drug_name: str) -> FDALabelInfo | None:
        """Get cached FDA label info for a drug.

        Args:
            drug_name: Drug name (generic or brand)

        Returns:
            FDALabelInfo if cached, None otherwise
        """
        self._ensure_loaded()

        normalized = drug_name.lower().strip()

        cached_data = None
        # Try exact match
        if normalized in self._cache:
            cached_data = self._cache[normalized]
        else:
            # Try without common suffixes
            for suffix in [" hydrochloride", " hcl", " mesylate", " tosylate", " maleate"]:
                if normalized.endswith(suffix):
                    stripped = normalized[:-len(suffix)]
                    if stripped in self._cache:
                        cached_data = self._cache[stripped]
                        break

        if cached_data is None:
            return None

        # Create FDALabelInfo from cached data
        info = FDALabelInfo(**cached_data)

        # Generate PubMed URL if not present (for backwards compatibility with old cache)
        if not info.pubmed_url:
            from urllib.parse import quote
            search_term = info.generic_name or info.drug_name
            info.pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/?term={quote(search_term)}+cancer"

        return info

    def set(self, drug_name: str, info: FDALabelInfo) -> None:
        """Cache FDA label info for a drug."""
        self._ensure_loaded()
        normalized = drug_name.lower().strip()
        self._cache[normalized] = info.model_dump()

    def save(self) -> None:
        """Save cache to disk."""
        self.cache_path.parent.mkdir(parents=True, exist_ok=True)

        data = {
            "last_updated": datetime.now().isoformat(),
            "drugs": self._cache
        }

        with open(self.cache_path, "w") as f:
            json.dump(data, f, indent=2, default=str)

        logger.info(f"Saved FDA label cache: {len(self._cache)} drugs to {self.cache_path}")

    @property
    def last_updated(self) -> str | None:
        """Get the last update timestamp."""
        self._ensure_loaded()
        return self._last_updated

    def __len__(self) -> int:
        self._ensure_loaded()
        return len(self._cache)


# Module-level singleton for FDA label cache
_fda_label_cache: FDALabelCache | None = None


def get_fda_label_cache() -> FDALabelCache:
    """Get the singleton FDALabelCache instance."""
    global _fda_label_cache
    if _fda_label_cache is None:
        _fda_label_cache = FDALabelCache()
    return _fda_label_cache


def get_cached_fda_label(drug_name: str) -> FDALabelInfo | None:
    """Get FDA label info from cache (no API call)."""
    return get_fda_label_cache().get(drug_name)


def get_fda_label_with_fallback(drug_name: str, fetch_if_missing: bool = False) -> FDALabelInfo | None:
    """Get FDA label info from cache, optionally fetching if missing.

    Args:
        drug_name: Drug name (generic or brand)
        fetch_if_missing: If True, fetch from API if not in cache (adds latency)

    Returns:
        FDALabelInfo if found/fetched, None otherwise
    """
    cached = get_fda_label_cache().get(drug_name)
    if cached:
        return cached

    if not fetch_if_missing:
        return None

    # Fetch on-demand
    import asyncio
    try:
        loop = asyncio.get_running_loop()
        # Can't use asyncio.run() in async context
        return None  # Skip fetch if already in async context
    except RuntimeError:
        # No running loop, safe to fetch
        try:
            result = asyncio.run(get_fda_label_info(drug_name))
            if result.api_success:
                get_fda_label_cache().set(drug_name, result)
                get_fda_label_cache().save()
            return result
        except Exception as e:
            logger.warning(f"Failed to fetch FDA label for {drug_name}: {e}")
            return None


async def refresh_fda_label_cache(
    drugs: list[str] | None = None,
    save: bool = True
) -> dict[str, FDALabelInfo]:
    """Refresh FDA label cache by fetching from OpenFDA API.

    This is intended to be run as a nightly batch job.

    Args:
        drugs: List of drug names to fetch. If None, fetches all drugs in ONCOLOGY_DRUG_MAPPINGS.
        save: Whether to save the cache to disk after fetching.

    Returns:
        Dictionary of drug_name -> FDALabelInfo
    """
    import asyncio

    cache = get_fda_label_cache()
    client = get_openfda_client()

    # If no drugs specified, use all from curated mapping
    if drugs is None:
        # Get unique generic names from the mapping
        seen_generics = set()
        drugs = []
        for drug_name, info in ONCOLOGY_DRUG_MAPPINGS.items():
            generic = info.get("generic", drug_name)
            if generic.lower() not in seen_generics:
                seen_generics.add(generic.lower())
                drugs.append(generic)

    logger.info(f"Refreshing FDA label cache for {len(drugs)} drugs...")

    results = {}

    # Fetch in batches to avoid overwhelming the API
    batch_size = 5
    for i in range(0, len(drugs), batch_size):
        batch = drugs[i:i + batch_size]
        tasks = [client.get_label_info(drug) for drug in batch]
        batch_results = await asyncio.gather(*tasks, return_exceptions=True)

        for drug, result in zip(batch, batch_results):
            if isinstance(result, Exception):
                logger.warning(f"Failed to fetch FDA label for {drug}: {result}")
            else:
                cache.set(drug, result)
                results[drug] = result
                logger.debug(f"Cached FDA label for {drug}")

        # Small delay between batches to be nice to the API
        if i + batch_size < len(drugs):
            await asyncio.sleep(0.5)

    if save:
        cache.save()

    logger.info(f"FDA label cache refresh complete: {len(results)} drugs cached")
    return results


async def ensure_fda_labels_cached(
    drugs: list[str],
    save: bool = True
) -> dict[str, FDALabelInfo]:
    """Ensure FDA labels are cached for a list of drugs.

    Only fetches drugs that are not already in the cache.
    Use this when you need FDA label data but haven't run the nightly refresh.

    Args:
        drugs: List of drug names to ensure are cached.
        save: Whether to save the cache after fetching missing drugs.

    Returns:
        Dictionary of drug_name -> FDALabelInfo for all requested drugs.
    """
    import asyncio

    cache = get_fda_label_cache()
    client = get_openfda_client()

    results = {}
    missing = []

    # Check which drugs are already cached
    for drug in drugs:
        cached = cache.get(drug)
        if cached:
            results[drug] = cached
        else:
            missing.append(drug)

    if not missing:
        logger.debug(f"All {len(drugs)} drugs already cached")
        return results

    logger.info(f"Fetching {len(missing)} missing drugs (out of {len(drugs)} requested)...")

    # Fetch missing drugs in batches
    batch_size = 5
    for i in range(0, len(missing), batch_size):
        batch = missing[i:i + batch_size]
        tasks = [client.get_label_info(drug) for drug in batch]
        batch_results = await asyncio.gather(*tasks, return_exceptions=True)

        for drug, result in zip(batch, batch_results):
            if isinstance(result, Exception):
                logger.warning(f"Failed to fetch FDA label for {drug}: {result}")
            else:
                cache.set(drug, result)
                results[drug] = result

        # Small delay between batches
        if i + batch_size < len(missing):
            await asyncio.sleep(0.3)

    if save and missing:
        cache.save()

    logger.info(f"FDA label cache updated: {len(results)} drugs now cached")
    return results


def ensure_fda_labels_cached_sync(drugs: list[str], save: bool = True) -> dict[str, FDALabelInfo]:
    """Synchronous wrapper for ensure_fda_labels_cached.

    Use this when you need to fetch missing drugs from synchronous code.
    """
    import asyncio

    try:
        loop = asyncio.get_running_loop()
        # If we're in an async context, can't use run()
        import concurrent.futures
        with concurrent.futures.ThreadPoolExecutor() as executor:
            future = executor.submit(asyncio.run, ensure_fda_labels_cached(drugs, save))
            return future.result()
    except RuntimeError:
        # No running loop, safe to use asyncio.run()
        return asyncio.run(ensure_fda_labels_cached(drugs, save))
