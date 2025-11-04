
.. |br| raw:: html

     <br>

############
Installation
############

Setting up ``pandeia``
======================

``pandeia`` is the ``python`` engine behind the JWST exposure time calculator
(ETC) (`Pontoppidan+ 2016 <https://ui.adsabs.harvard.edu/abs/2016SPIE.9910E..16P>`_).
``jayrock`` uses it to simulate JWST observations locally on your machine, avoiding
the need to use the online ETC.

``pandeia`` requires different reference datasets. Detailed install
instructions are `here
<https://outerspace.stsci.edu/spaces/PEN/pages/77530136/Pandeia+Engine+Installation>`_.
In short, you have to download three datasets and tell your system where to
find them by setting two system variables in your ``.bashrc`` or ``.zshrc`` (← Mac users) files.

1. JWST reference data (4GB) → `Section 3.1 <https://outerspace.stsci.edu/spaces/PEN/pages/77530136/Pandeia+Engine+Installation#PandeiaEngineInstallation-RequiredData>`_ of the instructions

.. code-block:: bash

   export pandeia_refdata=/path/to/jwst_refdata/  # in your .bashrc or .zshrc file

2. Synphot datasets 5 (Phoenix models, 1.7GB) and 7 (JWST ETC spectra, 8.5MB) → `Section 3.2 <https://outerspace.stsci.edu/spaces/PEN/pages/77530136/Pandeia+Engine+Installation#PandeiaEngineInstallation-RecommendedData>`_

.. code-block:: bash

   export PYSYN_CDBS=/path/to/synphot_datasets/  # in your .bashrc or .zshrc file

.. dropdown:: Example configuration

    Getting the directory structure of the synphot datasets right can be tricky. For reference, this is what it looks like on my system:

    .. code-block:: bash

       $ export pandeia_refdata=$HOME/astro/data/pandeia  # this line should be in your .zshrc or .bashrc file
       $ export PYSYN_CDBS=$HOME/astro/data/synphot  # this line should be in your .zshrc or .bashrc file
       $ ls $PYSYN_CDB
       comp  grid  grp
       $ ls $PYSYN_CDBS/comp
       acs  cos  fgs  foc  fos  hrs  hsp  nicmos  nonhst  ota  stis  wfc3  wfpc  wfpc2
       $ ls $PYSYN_CDBS/grid
       atmo2020  brown19  comp_qso  novae  phoenix  pne  solsys  stellar_pop  swire
       $ ls $PYSYN_CDBS/grp
       redcat



Installing ``jayrock``
======================

``jayrock`` requires ``python`` 3.11 or higher. Install it via

.. prompt:: bash
   :prompts: $

   pip install jayrock  # lines starting with $ are terminal commands

This will also install all required dependencies, including common ones like ``numpy``, ``pandas``,
``astropy``, ``astroquery``, and less common ones like ``rocks``, the `client <https://rocks.readthedocs.io>`_ to the `SsODNet database <https://ssp.imcce.fr/webservices/ssodnet/>`_ of minor body data, and ``jwst_gtvt``, the `JWST General Target Visibility Tool <https://jwst-docs.stsci.edu/jwst-other-tools/jwst-target-visibility-tools/jwst-general-target-visibility-tool-help>`_.

After installing ``jayrock``, you can check that it installed correctly by running

.. prompt:: python
   :prompts: >>>

   import jayrock  # lines starting with >>> are python commands

On your first run of a ``jayrock`` script, ``rocks`` might ask to `download data to a cache directory <https://rocks.readthedocs.io/en/latest/ssodnet.html#data-stored-on-your-machine>`_.
Confirm by typing ``y`` and pressing enter.
