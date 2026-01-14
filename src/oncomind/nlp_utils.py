# nlp_utils.py
import spacy
from negspacy.negation import Negex

# Load once when module is imported
nlp = spacy.load("en_core_web_sm")
nlp.add_pipe("negex")

diminishers = ["less", "lower", "reduced", "worse", "limited", "poor", "weak"]

def is_positive(statement: str, keywords: list[str]) -> bool:
    doc = nlp(statement)

    for token in doc:
        if any(token.lemma_.lower().startswith(kw.lower()) for kw in keywords):
            for child in token.head.children:
                if child.dep_ == "neg":
                    return False

            # Check if diminisher modifies this keyword
            for child in token.children:
                if child.text.lower() in diminishers:
                    return False

            return True

    return False

def is_positive_debug(statement: str, keywords: list[str]) -> bool:
    doc = nlp(statement)

    for token in doc:
        if any(token.lemma_.lower().startswith(kw.lower()) for kw in keywords):
            if lemma.startswith(("in", "un", "non", "im")) and lemma != kw:
                return False
            print(f"MATCHED: '{token.text}' (lemma: {token.lemma_})")
            print(f"  HEAD: '{token.head.text}'")
            for child in token.head.children:
                print(f"    CHILD: '{child.text}' dep={child.dep_}")
                if child.dep_ == "neg":
                    print("  → NEGATION FOUND, returning False")
                    return False
            return True

    print("No keyword matched")
    return False