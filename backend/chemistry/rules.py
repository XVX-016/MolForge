
from typing import List, Optional
from ..schemas.rules import ChemistryRule, RuleOperation, RuleMetadata, PropertyImpact, RuleSet

# Registry of all deterministic chemistry rules
# These are the ONLY rules the AI is allowed to select from.
OPTIMIZATION_RULES = [
    ChemistryRule(
        id="opt-fluorinate",
        title="Aromatic Fluorination",
        description="Replace a C-H with C-F to adjust metabolic stability and modulate LogP.",
        operation=RuleOperation.REPLACE,
        smarts_pattern="[c;H1]",
        replacement_smiles="F",
        impact=PropertyImpact(logp_delta=0.1, tpsa_delta=0.0, mw_delta=18.0),
        metadata=RuleMetadata(
            rationale="Fluorine is a common bioisostere that can increase metabolic stability without massive volume increase.",
            typical_use="When metabolic stability is low or aromatic rings need electronic tuning.",
            tags=["metabolic-stability", "logp", "bioisostere"]
        )
    ),
    ChemistryRule(
        id="opt-pyridine",
        title="Phenyl to Pyridine",
        description="Replace a phenyl ring with a pyridine ring to reduce LogP and increase solubility.",
        operation=RuleOperation.REPLACE,
        smarts_pattern="c1ccccc1",
        replacement_smiles="c1ccncc1",
        impact=PropertyImpact(logp_delta=-0.8, tpsa_delta=12.9, mw_delta=1.0),
        metadata=RuleMetadata(
            rationale="Inserting a nitrogen atom reduces lipophilicity and increases polar character.",
            typical_use="When LogP is too high and compound has poor solubility.",
            tags=["solubility", "logp", "bioisostere"]
        )
    ),
    ChemistryRule(
        id="opt-hydroxylate",
        title="Hydroxylation",
        description="Add a hydroxyl (-OH) group to increase polar surface area and reduce lipophilicity.",
        operation=RuleOperation.ADD,
        smarts_pattern="[C,c]",
        replacement_smiles="O",
        impact=PropertyImpact(logp_delta=-0.7, tpsa_delta=20.2, mw_delta=16.0),
        metadata=RuleMetadata(
            rationale="Directly increases solubility and H-bonding potential.",
            typical_use="When TPSA is low and solubility is a concern.",
            tags=["solubility", "tpsa", "hydrogen-bond"]
        )
    ),
    ChemistryRule(
        id="opt-methylate",
        title="Methylation",
        description="Add a methyl group to fill a hydrophobic pocket or block metabolic sites.",
        operation=RuleOperation.ADD,
        smarts_pattern="[C,c,N,O]",
        replacement_smiles="C",
        impact=PropertyImpact(logp_delta=0.5, tpsa_delta=0.0, mw_delta=14.0),
        metadata=RuleMetadata(
            rationale="Methyl groups can increase potency by filling small hydrophobic pockets (magic methyl effect).",
            typical_use="When potency needs a slight boost and LogP capacity remains.",
            tags=["potency", "hydrophobicity", "magic-methyl"]
        )
    )
]

def get_rule_by_id(rule_id: str) -> Optional[ChemistryRule]:
    for rule in OPTIMIZATION_RULES:
        if rule.id == rule_id:
            return rule
    return None

def get_applicable_rules(mol_smiles: str) -> List[ChemistryRule]:
    # Placeholder for smarter heuristic-based filtering
    # In a real implementation, this would check SMARTS matches and property profiles
    return OPTIMIZATION_RULES
