import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
from typing import Dict, Tuple
import seaborn as sns

class NanoplateVisualizer:
    """Create nanoplate-style visualizations for ddPCR results."""
    
    COLORS = {
        'negative': '#000000',
        'fam_only': '#FF00FF',
        'vic_only': '#00FF00',
        'double_positive': '#FFFFFF',
        'grid_line': '#333333',
        'border': '#FFFFFF'
    }
    
    def __init__(self, n_droplets: int = 20000):
        """Initialize nanoplate visualizer with droplet count."""
        self.n_droplets = n_droplets
        self.grid_shape = self._calculate_grid_dimensions()
    
    def _calculate_grid_dimensions(self) -> Tuple[int, int]:
        """Calculate optimal grid dimensions for display (rows, cols)."""
        aspect_ratio = 1.3
        rows = int(np.sqrt(self.n_droplets / aspect_ratio))
        cols = int(np.ceil(self.n_droplets / rows))
        return rows, cols
    
    def create_partition_array(self, data_dict: Dict) -> np.ndarray:
        """Create 2D array representing partition status (0=negative, 1=FAM, 2=VIC, 3=both)."""
        n_droplets = data_dict['n_droplets']
        fam_positive = data_dict['fam_positive']
        vic_positive = data_dict['vic_positive']
        partition_status = np.zeros(n_droplets, dtype=int)
        partition_status[fam_positive & ~vic_positive] = 1
        partition_status[~fam_positive & vic_positive] = 2
        partition_status[fam_positive & vic_positive] = 3
        rows, cols = self.grid_shape
        total_cells = rows * cols
        if n_droplets < total_cells:
            partition_status = np.pad(partition_status, (0, total_cells - n_droplets), constant_values=-1)
        elif n_droplets > total_cells:
            partition_status = partition_status[:total_cells]
        return partition_status.reshape(rows, cols)
    
    # def plot_nanoplate(self, data_dict: Dict, title: str = None, 
    #                   figsize: Tuple[int, int] = (12, 8)) -> Tuple:
    #     """Generate nanoplate visualization with partition grid."""
    #     partition_array = self.create_partition_array(data_dict)
    #     fig, ax = plt.subplots(figsize=figsize, facecolor='black')
    #     ax.set_facecolor('black')
    #     colors_list = ['#000000', self.COLORS['negative'], self.COLORS['fam_only'], 
    #                   self.COLORS['vic_only'], self.COLORS['double_positive']]
    #     cmap = ListedColormap(colors_list)
    #     im = ax.imshow(partition_array, cmap=cmap, aspect='equal',
    #                   vmin=-1, vmax=3, interpolation='nearest')
    #     rows, cols = partition_array.shape
    #     if rows < 100:
    #         for i in range(rows + 1):
    #             ax.axhline(i - 0.5, color=self.COLORS['grid_line'], linewidth=0.5, alpha=0.3)
    #         for j in range(cols + 1):
    #             ax.axvline(j - 0.5, color=self.COLORS['grid_line'], linewidth=0.5, alpha=0.3)
    #     ax.set_xticks([])
    #     ax.set_yticks([])
    #     for spine in ax.spines.values():
    #         spine.set_edgecolor(self.COLORS['border'])
    #         spine.set_linewidth(2)
    #     legend_elements = [
    #         mpatches.Patch(facecolor=self.COLORS['negative'], edgecolor='white', linewidth=0.5,
    #                       label=f"Negative ({data_dict['n_negative']:,})"),
    #         mpatches.Patch(facecolor=self.COLORS['fam_only'], edgecolor='white', linewidth=0.5,
    #                       label=f"Mutant Only ({data_dict['n_fam_only']:,})"),
    #         mpatches.Patch(facecolor=self.COLORS['vic_only'], edgecolor='white', linewidth=0.5,
    #                       label=f"Wildtype Only ({data_dict['n_vic_only']:,})"),
    #         mpatches.Patch(facecolor=self.COLORS['double_positive'], edgecolor='black', linewidth=0.5,
    #                       label=f"Both Targets ({data_dict['n_double_positive']:,})")
    #     ]
    #     legend = ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1), 
    #              frameon=True, fontsize=10, title='Partition Status', title_fontsize=11, 
    #              facecolor='#1a1a1a', edgecolor='white')
    #     plt.setp(legend.get_texts(), color='white')
    #     plt.setp(legend.get_title(), color='white')
    #     if title:
    #         title_text = title
    #     else:
    #         vaf_measured = data_dict.get('vaf_measured', 0)
    #         title_text = f"ddPCR Nanoplate - VAF: {vaf_measured:.3f}%"
    #     ax.set_title(title_text, fontsize=14, fontweight='bold', pad=20, color='white')
    #     total_positive = (data_dict['n_fam_only'] + data_dict['n_vic_only'] + 
    #                      data_dict['n_double_positive'])
    #     summary_text = (f"Total Partitions: {data_dict['n_droplets']:,}\n"
    #                    f"Positive: {total_positive:,} ({100*total_positive/data_dict['n_droplets']:.1f}%)\n"
    #                    f"Negative: {data_dict['n_negative']:,}")
    #     ax.text(0.02, 0.98, summary_text, transform=ax.transAxes, fontsize=9,
    #            verticalalignment='top', color='white',
    #            bbox=dict(boxstyle='round', facecolor='#1a1a1a', alpha=0.9, edgecolor='white'))
    #     plt.tight_layout()
    #     return fig, ax
    def plot_nanoplate(self, data_dict: Dict, title: str = None, 
                    figsize: Tuple[int, int] = (12, 8)) -> Tuple:
        """Generate nanoplate visualization with partition grid."""
        partition_array = self.create_partition_array(data_dict)
        fig, ax = plt.subplots(figsize=figsize, facecolor='black')
        ax.set_facecolor('black')
        colors_list = ['#000000', self.COLORS['negative'], self.COLORS['fam_only'], 
                    self.COLORS['vic_only'], self.COLORS['double_positive']]
        cmap = ListedColormap(colors_list)
        im = ax.imshow(partition_array, cmap=cmap, aspect='equal',
                    vmin=-1, vmax=3, interpolation='nearest')
        rows, cols = partition_array.shape
        if rows < 100:
            for i in range(rows + 1):
                ax.axhline(i - 0.5, color=self.COLORS['grid_line'], linewidth=0.5, alpha=0.3)
            for j in range(cols + 1):
                ax.axvline(j - 0.5, color=self.COLORS['grid_line'], linewidth=0.5, alpha=0.3)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_edgecolor(self.COLORS['border'])
            spine.set_linewidth(2)
        
        # Use n_double_negative if n_negative doesn't exist
        n_negative = data_dict.get('n_negative', data_dict.get('n_double_negative', 0))
        
        legend_elements = [
            mpatches.Patch(facecolor=self.COLORS['negative'], edgecolor='white', linewidth=0.5,
                        label=f"Negative ({n_negative:,})"),
            mpatches.Patch(facecolor=self.COLORS['fam_only'], edgecolor='white', linewidth=0.5,
                        label=f"Mutant Only ({data_dict['n_fam_only']:,})"),
            mpatches.Patch(facecolor=self.COLORS['vic_only'], edgecolor='white', linewidth=0.5,
                        label=f"Wildtype Only ({data_dict['n_vic_only']:,})"),
            mpatches.Patch(facecolor=self.COLORS['double_positive'], edgecolor='black', linewidth=0.5,
                        label=f"Both Targets ({data_dict['n_double_positive']:,})")
        ]
        legend = ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1), 
                frameon=True, fontsize=10, title='Partition Status', title_fontsize=11, 
                facecolor='#1a1a1a', edgecolor='white')
        plt.setp(legend.get_texts(), color='white')
        plt.setp(legend.get_title(), color='white')
        if title:
            title_text = title
        else:
            vaf_measured = data_dict.get('vaf_measured', 0)
            title_text = f"ddPCR Nanoplate - VAF: {vaf_measured:.3f}%"
        ax.set_title(title_text, fontsize=14, fontweight='bold', pad=20, color='white')
        total_positive = (data_dict['n_fam_only'] + data_dict['n_vic_only'] + 
                        data_dict['n_double_positive'])
        summary_text = (f"Total Partitions: {data_dict['n_droplets']:,}\n"
                    f"Positive: {total_positive:,} ({100*total_positive/data_dict['n_droplets']:.1f}%)\n"
                    f"Negative: {n_negative:,}")
        ax.text(0.02, 0.98, summary_text, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', color='white',
            bbox=dict(boxstyle='round', facecolor='#1a1a1a', alpha=0.9, edgecolor='white'))
        plt.tight_layout()
        return fig, ax

    def plot_nanoplate_zoom(self, data_dict: Dict, 
                           zoom_region: Tuple[int, int, int, int] = None,
                           figsize: Tuple[int, int] = (10, 10)) -> Tuple:
        """Create zoomed-in view of nanoplate region showing individual partitions."""
        partition_array = self.create_partition_array(data_dict)
        if zoom_region is None:
            rows, cols = partition_array.shape
            zoom_size = min(50, rows//2)
            row_start = (rows - zoom_size) // 2
            col_start = (cols - zoom_size) // 2
            zoom_region = (row_start, row_start + zoom_size, col_start, col_start + zoom_size)
        row_start, row_end, col_start, col_end = zoom_region
        zoomed_array = partition_array[row_start:row_end, col_start:col_end]
        fig, ax = plt.subplots(figsize=figsize, facecolor='black')
        ax.set_facecolor('black')
        colors_list = ['#000000', self.COLORS['negative'], self.COLORS['fam_only'],
                      self.COLORS['vic_only'], self.COLORS['double_positive']]
        cmap = ListedColormap(colors_list)
        im = ax.imshow(zoomed_array, cmap=cmap, aspect='equal', vmin=-1, vmax=3, interpolation='nearest')
        zoom_rows, zoom_cols = zoomed_array.shape
        for i in range(zoom_rows + 1):
            ax.axhline(i - 0.5, color=self.COLORS['grid_line'], linewidth=1)
        for j in range(zoom_cols + 1):
            ax.axvline(j - 0.5, color=self.COLORS['grid_line'], linewidth=1)
        if zoom_rows <= 20 and zoom_cols <= 20:
            label_map = {0: 'N', 1: 'M', 2: 'W', 3: 'B'}
            text_colors = {0: 'white', 1: 'black', 2: 'black', 3: 'black'}
            for i in range(zoom_rows):
                for j in range(zoom_cols):
                    value = zoomed_array[i, j]
                    if value >= 0:
                        ax.text(j, i, label_map[value], ha='center', va='center',
                               fontsize=8, fontweight='bold', color=text_colors[value])
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_edgecolor(self.COLORS['border'])
            spine.set_linewidth(2)
        vaf_measured = data_dict.get('vaf_measured', 0)
        ax.set_title(f'Zoomed Nanoplate Region - VAF: {vaf_measured:.3f}%',
                    fontsize=14, fontweight='bold', pad=15, color='white')
        legend_elements = [
            mpatches.Patch(facecolor=self.COLORS['negative'], edgecolor='white', label='N: Negative'),
            mpatches.Patch(facecolor=self.COLORS['fam_only'], edgecolor='white', label='M: Mutant'),
            mpatches.Patch(facecolor=self.COLORS['vic_only'], edgecolor='white', label='W: Wildtype'),
            mpatches.Patch(facecolor=self.COLORS['double_positive'], edgecolor='black', label='B: Both')
        ]
        legend = ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1),
                 frameon=True, facecolor='#1a1a1a', edgecolor='white')
        plt.setp(legend.get_texts(), color='white')
        plt.tight_layout()
        return fig, ax
    
    def plot_multi_sample_comparison(self, results_dict: Dict,
                                    figsize: Tuple[int, int] = (16, 10)) -> Tuple:
        """Compare multiple samples in grid layout."""
        n_samples = len(results_dict)
        n_cols = min(3, n_samples)
        n_rows = int(np.ceil(n_samples / n_cols))
        fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, facecolor='black')
        fig.patch.set_facecolor('black')
        if n_samples == 1:
            axes = np.array([axes])
        axes = axes.flatten()
        colors_list = ['#000000', self.COLORS['negative'], self.COLORS['fam_only'],
                      self.COLORS['vic_only'], self.COLORS['double_positive']]
        cmap = ListedColormap(colors_list)
        for idx, (sample_name, data_dict) in enumerate(results_dict.items()):
            ax = axes[idx]
            ax.set_facecolor('black')
            partition_array = self.create_partition_array(data_dict)
            ax.imshow(partition_array, cmap=cmap, aspect='equal', vmin=-1, vmax=3, interpolation='nearest')
            ax.set_xticks([])
            ax.set_yticks([])
            vaf = data_dict.get('vaf_measured', 0)
            ax.set_title(f'{sample_name}\nVAF: {vaf:.3f}%', fontsize=11, fontweight='bold', color='white')
            for spine in ax.spines.values():
                spine.set_edgecolor(self.COLORS['border'])
                spine.set_linewidth(1.5)
        for idx in range(n_samples, len(axes)):
            axes[idx].axis('off')
            axes[idx].set_facecolor('black')
        legend_elements = [
            mpatches.Patch(facecolor=self.COLORS['negative'], edgecolor='white', label='Negative'),
            mpatches.Patch(facecolor=self.COLORS['fam_only'], edgecolor='white', label='Mutant'),
            mpatches.Patch(facecolor=self.COLORS['vic_only'], edgecolor='white', label='Wildtype'),
            mpatches.Patch(facecolor=self.COLORS['double_positive'], edgecolor='black', label='Both')
        ]
        legend = fig.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, 0.98),
                  ncol=4, frameon=True, fontsize=11, facecolor='#1a1a1a', edgecolor='white')
        plt.setp(legend.get_texts(), color='white')
        plt.suptitle('ddPCR Nanoplate Comparison', fontsize=16, fontweight='bold', y=0.995, color='white')
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        return fig, axes


# def create_publication_figure(data_dict: Dict, output_path: str = None) -> Tuple:
#     """Create publication-quality figure combining nanoplate and amplitude plots."""
#     visualizer = NanoplateVisualizer(n_droplets=data_dict['n_droplets'])
#     fig = plt.figure(figsize=(16, 10), facecolor='black')
#     fig.patch.set_facecolor('black')
#     gs = fig.add_gridspec(2, 2, height_ratios=[1.2, 1], width_ratios=[1, 1], hspace=0.3, wspace=0.3)
#     ax_nano = fig.add_subplot(gs[0, :])
#     ax_nano.set_facecolor('black')
#     partition_array = visualizer.create_partition_array(data_dict)
#     colors_list = ['#000000', visualizer.COLORS['negative'], visualizer.COLORS['fam_only'],
#                   visualizer.COLORS['vic_only'], visualizer.COLORS['double_positive']]
#     cmap = ListedColormap(colors_list)
#     ax_nano.imshow(partition_array, cmap=cmap, aspect='equal', vmin=-1, vmax=3, interpolation='nearest')
#     ax_nano.set_xticks([])
#     ax_nano.set_yticks([])
#     for spine in ax_nano.spines.values():
#         spine.set_edgecolor(visualizer.COLORS['border'])
#         spine.set_linewidth(2)
#     vaf_measured = data_dict.get('vaf_measured', 0)
#     ax_nano.set_title(f'ddPCR Nanoplate Visualization - VAF: {vaf_measured:.3f}%',
#                      fontsize=14, fontweight='bold', pad=15, color='white')
#     legend_elements = [
#         mpatches.Patch(facecolor=visualizer.COLORS['negative'], edgecolor='white', linewidth=0.5,
#                       label=f"Negative ({data_dict['n_negative']:,})"),
#         mpatches.Patch(facecolor=visualizer.COLORS['fam_only'], edgecolor='white', linewidth=0.5,
#                       label=f"Mutant ({data_dict['n_fam_only']:,})"),
#         mpatches.Patch(facecolor=visualizer.COLORS['vic_only'], edgecolor='white', linewidth=0.5,
#                       label=f"Wildtype ({data_dict['n_vic_only']:,})"),
#         mpatches.Patch(facecolor=visualizer.COLORS['double_positive'], edgecolor='black', linewidth=0.5,
#                       label=f"Both ({data_dict['n_double_positive']:,})")
#     ]
#     legend = ax_nano.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1.02, 0.5),
#                   frameon=True, fontsize=10, facecolor='#1a1a1a', edgecolor='white')
#     plt.setp(legend.get_texts(), color='white')
#     ax_fam = fig.add_subplot(gs[1, 0])
#     ax_fam.set_facecolor('black')
#     fluorescence_fam = data_dict['fam_fluorescence']
#     positive_mask_fam = data_dict['fam_positive']
#     threshold_fam = 3500
#     event_numbers = np.arange(len(fluorescence_fam))
#     ax_fam.scatter(event_numbers[~positive_mask_fam], fluorescence_fam[~positive_mask_fam],
#                   c='#404040', s=2, alpha=0.3, label='Negative')
#     ax_fam.scatter(event_numbers[positive_mask_fam], fluorescence_fam[positive_mask_fam],
#                   c='#FF00FF', s=2, alpha=0.6, label='Positive')
#     ax_fam.axhline(threshold_fam, color='red', linestyle='--', linewidth=2, label='Threshold')
#     ax_fam.set_xlabel('Event Number', fontweight='bold', color='white')
#     ax_fam.set_ylabel('FAM Amplitude', fontweight='bold', color='white')
#     ax_fam.set_title('FAM Channel (Mutant)', fontweight='bold', color='white')
#     legend_fam = ax_fam.legend(loc='upper right', fontsize=8, facecolor='#1a1a1a', edgecolor='white')
#     plt.setp(legend_fam.get_texts(), color='white')
#     ax_fam.grid(True, alpha=0.3, color='#404040')
#     ax_fam.tick_params(colors='white')
#     for spine in ax_fam.spines.values():
#         spine.set_edgecolor('white')
#     ax_vic = fig.add_subplot(gs[1, 1])
#     ax_vic.set_facecolor('black')
#     fluorescence_vic = data_dict['vic_fluorescence']
#     positive_mask_vic = data_dict['vic_positive']
#     threshold_vic = 3500
#     ax_vic.scatter(event_numbers[~positive_mask_vic], fluorescence_vic[~positive_mask_vic],
#                   c='#404040', s=2, alpha=0.3, label='Negative')
#     ax_vic.scatter(event_numbers[positive_mask_vic], fluorescence_vic[positive_mask_vic],
#                   c='#00FF00', s=2, alpha=0.6, label='Positive')
#     ax_vic.axhline(threshold_vic, color='red', linestyle='--', linewidth=2, label='Threshold')
#     ax_vic.set_xlabel('Event Number', fontweight='bold', color='white')
#     ax_vic.set_ylabel('VIC Amplitude', fontweight='bold', color='white')
#     ax_vic.set_title('VIC Channel (Wildtype)', fontweight='bold', color='white')
#     legend_vic = ax_vic.legend(loc='upper right', fontsize=8, facecolor='#1a1a1a', edgecolor='white')
#     plt.setp(legend_vic.get_texts(), color='white')
#     ax_vic.grid(True, alpha=0.3, color='#404040')
#     ax_vic.tick_params(colors='white')
#     for spine in ax_vic.spines.values():
#         spine.set_edgecolor('white')
#     if output_path:
#         plt.savefig(output_path, dpi=300, bbox_inches='tight')
#         print(f"Figure saved: {output_path}")
#     return fig, (ax_nano, ax_fam, ax_vic)


def create_publication_figure(data_dict: Dict, output_path: str = None) -> Tuple:
    """Create publication-quality figure combining nanoplate and amplitude plots."""
    visualizer = NanoplateVisualizer(n_droplets=data_dict['n_droplets'])
    fig = plt.figure(figsize=(16, 10), facecolor='black')
    fig.patch.set_facecolor('black')
    gs = fig.add_gridspec(2, 2, height_ratios=[1.2, 1], width_ratios=[1, 1], hspace=0.3, wspace=0.3)
    ax_nano = fig.add_subplot(gs[0, :])
    ax_nano.set_facecolor('black')
    partition_array = visualizer.create_partition_array(data_dict)
    colors_list = ['#000000', visualizer.COLORS['negative'], visualizer.COLORS['fam_only'],
                  visualizer.COLORS['vic_only'], visualizer.COLORS['double_positive']]
    cmap = ListedColormap(colors_list)
    ax_nano.imshow(partition_array, cmap=cmap, aspect='equal', vmin=-1, vmax=3, interpolation='nearest')
    ax_nano.set_xticks([])
    ax_nano.set_yticks([])
    for spine in ax_nano.spines.values():
        spine.set_edgecolor(visualizer.COLORS['border'])
        spine.set_linewidth(2)
    vaf_measured = data_dict.get('vaf_measured', 0)
    ax_nano.set_title(f'ddPCR Nanoplate Visualization - VAF: {vaf_measured:.3f}%',
                     fontsize=14, fontweight='bold', pad=15, color='white')
    
    # Use n_double_negative if n_negative doesn't exist
    n_negative = data_dict.get('n_negative', data_dict.get('n_double_negative', 0))
    
    legend_elements = [
        mpatches.Patch(facecolor=visualizer.COLORS['negative'], edgecolor='white', linewidth=0.5,
                      label=f"Negative ({n_negative:,})"),
        mpatches.Patch(facecolor=visualizer.COLORS['fam_only'], edgecolor='white', linewidth=0.5,
                      label=f"Mutant ({data_dict['n_fam_only']:,})"),
        mpatches.Patch(facecolor=visualizer.COLORS['vic_only'], edgecolor='white', linewidth=0.5,
                      label=f"Wildtype ({data_dict['n_vic_only']:,})"),
        mpatches.Patch(facecolor=visualizer.COLORS['double_positive'], edgecolor='black', linewidth=0.5,
                      label=f"Both ({data_dict['n_double_positive']:,})")
    ]
    legend = ax_nano.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1.02, 0.5),
                  frameon=True, fontsize=10, facecolor='#1a1a1a', edgecolor='white')
    plt.setp(legend.get_texts(), color='white')
    
    # Add FAM channel amplitude plot
    ax_fam = fig.add_subplot(gs[1, 0])
    ax_fam.set_facecolor('black')
    fluorescence_fam = data_dict['fam_fluorescence']
    positive_mask_fam = data_dict['fam_positive']
    threshold_fam = 3500
    event_numbers = np.arange(len(fluorescence_fam))
    ax_fam.scatter(event_numbers[~positive_mask_fam], fluorescence_fam[~positive_mask_fam],
                  c='#404040', s=2, alpha=0.3, label='Negative')
    ax_fam.scatter(event_numbers[positive_mask_fam], fluorescence_fam[positive_mask_fam],
                  c='#FF00FF', s=2, alpha=0.6, label='Positive')
    ax_fam.axhline(threshold_fam, color='red', linestyle='--', linewidth=2, label='Threshold')
    ax_fam.set_xlabel('Event Number', fontweight='bold', color='white')
    ax_fam.set_ylabel('FAM Amplitude', fontweight='bold', color='white')
    ax_fam.set_title('FAM Channel (Mutant)', fontweight='bold', color='white')
    legend_fam = ax_fam.legend(loc='upper right', fontsize=8, facecolor='#1a1a1a', edgecolor='white')
    plt.setp(legend_fam.get_texts(), color='white')
    ax_fam.grid(True, alpha=0.3, color='#404040')
    ax_fam.tick_params(colors='white')
    for spine in ax_fam.spines.values():
        spine.set_edgecolor('white')
    
    # Add VIC channel amplitude plot
    ax_vic = fig.add_subplot(gs[1, 1])
    ax_vic.set_facecolor('black')
    fluorescence_vic = data_dict['vic_fluorescence']
    positive_mask_vic = data_dict['vic_positive']
    threshold_vic = 3500
    ax_vic.scatter(event_numbers[~positive_mask_vic], fluorescence_vic[~positive_mask_vic],
                  c='#404040', s=2, alpha=0.3, label='Negative')
    ax_vic.scatter(event_numbers[positive_mask_vic], fluorescence_vic[positive_mask_vic],
                  c='#00FF00', s=2, alpha=0.6, label='Positive')
    ax_vic.axhline(threshold_vic, color='red', linestyle='--', linewidth=2, label='Threshold')
    ax_vic.set_xlabel('Event Number', fontweight='bold', color='white')
    ax_vic.set_ylabel('VIC Amplitude', fontweight='bold', color='white')
    ax_vic.set_title('VIC Channel (Wildtype)', fontweight='bold', color='white')
    legend_vic = ax_vic.legend(loc='upper right', fontsize=8, facecolor='#1a1a1a', edgecolor='white')
    plt.setp(legend_vic.get_texts(), color='white')
    ax_vic.grid(True, alpha=0.3, color='#404040')
    ax_vic.tick_params(colors='white')
    for spine in ax_vic.spines.values():
        spine.set_edgecolor('white')
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved: {output_path}")
    
    return fig, (ax_nano, ax_fam, ax_vic)