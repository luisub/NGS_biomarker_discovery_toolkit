## !/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
by: Luis U. Aguilera
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pypcrtool.pcr import InSilicoPCR
from pathlib import Path

def PCRVisualization(primers, wildtype_fasta_path, mutant_fasta_path, mutation_name, output_directory, lane_label_list=None):
    """
    Simulate four PCR reactions (WT/Mut primers and WT/Mut templates) and draw a gel.
    Input: primers (dict), wildtype_fasta_path (Path), mutant_fasta_path (Path), mutation_name (str), output_directory (Path), lane_label_list (list[str] or None)
    Output: dict - mapping reaction label to product size in bp (or None if no product)
    """
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    reactions = [
        ("WT primers + WT template", primers["forward_wildtype"]["sequence"], primers["reverse_common"]["sequence"], wildtype_fasta_path),
        ("WT primers + mutant template", primers["forward_wildtype"]["sequence"], primers["reverse_common"]["sequence"], mutant_fasta_path),
        ("Mutant primers + WT template", primers["forward_mutant"]["sequence"], primers["reverse_common"]["sequence"], wildtype_fasta_path),
        ("Mutant primers + mutant template", primers["forward_mutant"]["sequence"], primers["reverse_common"]["sequence"], mutant_fasta_path),
    ]
    pcr_results = {}
    for index, (reaction_label, forward_sequence, reverse_sequence, template_path) in enumerate(reactions, start=1):
        print(f"{index}. {reaction_label}")
        insilico_pcr = InSilicoPCR(forward_primer=forward_sequence, reverse_primer=reverse_sequence, sequence_file=str(template_path))
        try:
            products = insilico_pcr.perform_pcr()
        except Exception as error:
            print(f"   Error during simulation: {error}")
            pcr_results[reaction_label] = None
            continue
        product_size = None
        if products and len(products) > 0:
            first_product = products[0]
            if isinstance(first_product, dict) and "amplicon" in first_product:
                amplicon = first_product["amplicon"]
                if isinstance(amplicon, str):
                    product_size = len(amplicon)
                elif isinstance(amplicon, (list, tuple)):
                    product_size = len(amplicon)
                elif isinstance(amplicon, (int, float)):
                    product_size = int(amplicon)
            elif isinstance(first_product, tuple):
                sequence_candidate = None
                for element in first_product:
                    if isinstance(element, str):
                        sequence_candidate = element
                        break
                if sequence_candidate is not None:
                    product_size = len(sequence_candidate)
                else:
                    for element in first_product:
                        if isinstance(element, (int, np.integer)) and element > 0:
                            product_size = int(element)
                            break
            elif isinstance(first_product, str):
                product_size = len(first_product)
            if product_size is not None:
                pcr_results[reaction_label] = product_size
                print(f"   Result: Product {product_size} bp")
            else:
                pcr_results[reaction_label] = None
                print("   Result: Product found, but could not determine size")
        else:
            pcr_results[reaction_label] = None
            print("   Result: No product")
    fig, ax = plt.subplots(figsize=(6, 10))
    ax.set_facecolor("#e0e0e0")
    fig.patch.set_facecolor("#e0e0e0")
    marker_lane_x = 0.7
    marker_sizes = [100, 150, 200, 300, 400, 500]
    for marker_size in marker_sizes:
        y_position = np.log10(marker_size) * 1.5
        band_rectangle = patches.Rectangle((marker_lane_x - 0.15, y_position - 0.01), 0.3, 0.02, facecolor="black", edgecolor="none", alpha=0.9)
        ax.add_patch(band_rectangle)
        ax.text(marker_lane_x - 0.3, y_position, f"{marker_size}", ha="right", va="center", color="black", fontsize=8)
    lane_x_positions = [2.0, 3.3, 4.6, 5.9]
    lanes_list = list(pcr_results.items())
    for lane_index, (lane_x, (reaction_label, product_size)) in enumerate(zip(lane_x_positions, lanes_list)):
        if lane_label_list is not None and lane_index < len(lane_label_list):
            lane_title = lane_label_list[lane_index]
        else:
            lane_title = reaction_label
        ax.text(lane_x, 4.3, lane_title, ha="center", va="bottom", color="black", fontsize=8, rotation=0, wrap=True)
        loading_well = patches.Rectangle((lane_x - 0.2, 4.0), 0.4, 0.05, facecolor="#b0b0b0", edgecolor="black", linewidth=0.8, alpha=1.0)
        ax.add_patch(loading_well)
        if product_size:
            y_position = np.log10(product_size) * 1.5
            product_band = patches.Rectangle((lane_x - 0.25, y_position - 0.01), 0.5, 0.02, facecolor="black", edgecolor="none", alpha=0.95)
            ax.add_patch(product_band)
    ax.set_xlim(0, 6.8)
    ax.set_ylim(0.5, 4.5)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    legend_elements = [patches.Patch(facecolor="black", alpha=0.9, label="DNA Ladder / PCR band")]
    legend = ax.legend(handles=legend_elements, loc="upper right", framealpha=0.9, facecolor="#f0f0f0", edgecolor="black", fontsize=8)
    plt.tight_layout()
    gel_image_path = output_directory / f"{mutation_name}_gel_electrophoresis.png"
    plt.savefig(gel_image_path, dpi=300, facecolor="#e0e0e0", bbox_inches="tight")
    print(f"Gel image saved: {gel_image_path}")
    plt.show()
    return pcr_results
