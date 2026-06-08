#!/usr/bin/env python3
"""Quick test: generate a reconstruction plot with the divergence panel."""
import sys, tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "src"))

from GeneFior.GeneFior_Reconstruct import GeneReconstructor, _analyze_divergence

def main():
    # Build synthetic data with ~12% divergence
    consensus = "ATCAATCTATCG" * 65   # 780 bp, differs at pos 4,5,10
    reference = "ATCGATCGATCG" * 65

    per_pos = [
        {
            "pos": i + 1,
            "depth": 18 + (i % 7),
            "consensus_base": consensus[i],
            "top_base_freq": 88.0 if consensus[i] != reference[i] else 97.0,
            "alt_base": reference[i] if consensus[i] != reference[i] else "",
            "alt_base_freq": 12.0 if consensus[i] != reference[i] else 0.0,
            "has_insertion": False,
            "ins_seq": "", "ins_depth": 0, "ins_freq": 0.0
        }
        for i in range(len(consensus))
    ]
    # Make a couple of positions "informative" (high alt freq)
    per_pos[99]["alt_base_freq"] = 44.0
    per_pos[99]["top_base_freq"] = 56.0

    informative_pos = [{"pos": 100}]
    variants = [
        {"pos": p["pos"], "type": "SNP",
         "ref_base": reference[p["pos"]-1],
         "sample_base": p["consensus_base"],
         "depth": p["depth"], "freq_pct": p["top_base_freq"]}
        for p in per_pos
        if p["consensus_base"] != reference[p["pos"]-1]
    ]

    val = {
        "grade": "B", "coverage_pct": 96.2, "mean_depth": 20.5,
        "read_support_pct": 91.4, "cv_depth_pct": 28.1,
        "gene_length": len(consensus)
    }

    div = _analyze_divergence(consensus, reference, per_pos, variants)
    print(f"Divergence: {div['divergence_pct']}%  ({div['divergence_class']})")
    print(f"SNPs: {div['snps']}")

    with tempfile.TemporaryDirectory() as d:
        td = Path(d)

        rec = GeneReconstructor.__new__(GeneReconstructor)
        rec.gene_name = "tet(Q)_1_L33696"
        rec.min_depth = 5

        rec._plot_reconstruction(
            tool_label="Bowtie2",
            tool_dir=td,
            gene_safe="tet_Q__1_L33696",
            per_pos=per_pos,
            variants=variants,
            informative_pos=informative_pos,
            haplotypes=[],
            mv_class="single",
            val=val,
            reference=reference,
            divergence=div,
        )

        out = td / "tet_Q__1_L33696_reconstruction_plot.html"
        assert out.exists(), "HTML not created"
        html = out.read_text()

        checks = {
            "divergence panel label": "Sample vs Reference divergence" in html,
            "divergence % in title":  f"{div['divergence_pct']:.1f}%" in html,
            "SNP colour present":     "#d73027" in html,
            "match colour present":   "#4dac26" in html,
            "biology note shown":     "EXPECTED BIOLOGY" in html,
            "divergence badge":       "divergence" in html,
        }

        all_ok = True
        for check, result in checks.items():
            sym = "✓" if result else "✗"
            print(f"  {sym} {check}")
            all_ok = all_ok and result

        print(f"\nHTML size: {len(html):,} bytes")
        print("PASS" if all_ok else "FAIL")
        return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())

