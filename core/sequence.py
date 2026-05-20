"""
CaloXSequence: four-step analysis unit bundling variable definition, event
selection, histogram booking, and plot generation under a single name.

Steps (all optional):
  define_vars(ctx)       -> rdf | None   Add derived RDF columns; return updated
                                         rdf or None to leave ctx.rdf unchanged.
  define_selection(ctx)  -> rdf | None   Return a filtered RDF view for this
                                         sequence only; does not affect other
                                         sequences.
  book_hists(ctx)        -> None          Book lazy ROOT histograms into
                                         ctx.hbook.  Called with ctx.rdf
                                         already set to the sequence-specific
                                         filtered view (if define_selection
                                         is provided).
  make_plots(ctx)        -> html_path(s) | None
                                         Generate HTML plots from the saved
                                         ROOT / JSON files.

Runner functions
----------------
  run_hist_phase(sequences, ctx)
      Step 1 — define_vars (sequential; ctx.rdf updated when a non-None rdf is
               returned)
      Step 2 — define_selection + book_hists (ctx.rdf is temporarily swapped to
               the per-sequence filtered view during booking; restored after)
      Step 3 — ROOT.RDF.RunGraphs (single call across all sequences)
      Step 4 — save callbacks (run after the graph is materialised)

  run_plot_phase(sequences, ctx)
      Calls make_plots for each sequence and prints a summary of HTML outputs.
"""

from dataclasses import dataclass
from typing import Callable, Optional
import ROOT


@dataclass
class CaloXSequence:
    name: str
    define_vars: Optional[Callable] = None       # (ctx) -> rdf | None
    define_selection: Optional[Callable] = None  # (ctx) -> rdf | None
    book_hists: Optional[Callable] = None        # (ctx) -> None  (registers into ctx.hbook)
    make_plots: Optional[Callable] = None        # (ctx) -> html_path(s) | None
    enabled_by_default: bool = True


def run_hist_phase(sequences, ctx):
    """Four-step hist phase: define_vars → define_selection → book_hists → RunGraphs → save."""
    # Step 1: define_vars — add derived columns sequentially
    for seq in sequences:
        if seq.define_vars is not None:
            result = seq.define_vars(ctx)
            if result is not None:
                ctx.rdf = result

    # Step 2: book_hists — each call registers into ctx.hbook
    for seq in sequences:
        if seq.book_hists is None:
            continue
        print(f"\033[94mBooking {seq.name}\033[0m")
        if seq.define_selection is not None:
            orig_rdf = ctx.rdf
            ctx.rdf = seq.define_selection(ctx)
            seq.book_hists(ctx)
            ctx.rdf = orig_rdf
        else:
            seq.book_hists(ctx)

    # Include cutflow count/sum proxies so they're evaluated in the same pass
    if hasattr(ctx, 'sel_mgr'):
        ctx.hbook.add(None, ctx.sel_mgr.all_cutflow_proxies)

    # Step 3: single RunGraphs call across all registered histograms
    lazy = ctx.hbook.lazy_objects
    if lazy:
        print(f"\n\033[94mTriggering {len(lazy)} lazy computations...\033[0m")
        ROOT.RDF.RunGraphs(lazy)

    # Step 4: save all ROOT files and run post-save callbacks
    print(f"\033[94mSaving all histograms...\033[0m")
    ctx.hbook.save_all()


def run_plot_phase(sequences, ctx):
    """Run make_plots for each sequence; collect and print HTML outputs."""
    outputs = {}
    for seq in sequences:
        if seq.make_plots is None:
            continue
        print(f"Generating {seq.name} plots...")
        outputs[seq.name] = seq.make_plots(ctx)

    print("\n" + "*" * 30)
    print("Plot Generation Summary:")
    for label, url in outputs.items():
        if isinstance(url, str):
            print(f"  {label}: {url}")
        elif isinstance(url, list):
            print(f"  {label}:")
            for u in url:
                print(f"    - {u}")

    return outputs
