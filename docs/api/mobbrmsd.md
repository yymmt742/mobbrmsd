# mobbRMSD

Python interface for molecular-oriented RMSD.

::: mobbrmsd.mobbrmsd.mobbrmsd
    options:
      members: true
      show_source: false
      show_docstring_parameters: true
      show_signature_annotations: false
      merge_init_into_class: true
      separate_signature: true
      show_signature: true
      inherited_members: false
      modernize_annotations: true
      members:
        - rmsd
        - run
        - batch_run
        - min_span_tree

::: mobbrmsd.mobbrmsd.mobbrmsd_result
    options:
      members: true
      show_source: false
      show_docstring_parameters: false
      show_signature_annotations: false
      merge_init_into_class: false
      separate_signature: true
      show_signature: true
      inherited_members: true
      modernize_annotations: true
      members:
        - autocorr
        - lowerbound
        - upperbound
        - lowerbound_as_rmsd
        - upperbound_as_rmsd
        - sd
        - rmsd
        - bounds
        - bounds_as_rmsd
        - n_eval
        - log_eval_ratio
        - eval_ratio
        - is_finished
        - restart
