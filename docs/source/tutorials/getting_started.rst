Getting Started
===============

This tutorial will guide you through the basic usage of protdata to load proteomics data into the AnnData format.

Basic Usage
-----------

protdata reads proteomics search engine output files into the `AnnData <https://anndata.readthedocs.io/>`_ format. All loaders return AnnData objects with a consistent structure:

- **X**: The main intensity matrix (samples × proteins)
- **obs**: Sample metadata
- **var**: Protein metadata
- **uns**: Unstructured metadata (search engine info, etc.)
- **layers**: Additional data matrices (e.g., spectral counts, raw intensities)

Loading MaxQuant Data
~~~~~~~~~~~~~~~~~~~~~

You can download an example proteinGroups `file here <https://zenodo.org/records/3774452/files/MaxQuant_Protein_Groups.tabular?download=1>`_

.. code-block:: python

   import protdata

   # Load MaxQuant data
   adata = protdata.io.read_maxquant("proteinGroups.txt")
   
   # Inspect the data
   print(adata)
   print(f"Samples: {adata.n_obs}")
   print(f"Proteins: {adata.n_vars}")
   print(f"Available layers: {list(adata.layers.keys())}")

Loading FragPipe Data
~~~~~~~~~~~~~~~~~~~~~

You can download an example FragPipe output `file here <https://github.com/Nesvilab/philosopher/blob/master/testdata/combined_protein.tsv?raw=true>`_

.. code-block:: python

   # Load FragPipe data
   adata = protdata.io.read_fragpipe("combined_protein.tsv")
   
   print(adata)

Loading DIA-NN Data
~~~~~~~~~~~~~~~~~~~

You can download an example DIA-NN report `file here <https://github.com/vdemichev/DiaNN/raw/master/diann-output-examples/report.pg_matrix.tsv>`_

.. code-block:: python

   # Load DIA-NN data
   adata = protdata.io.read_diann("report.pg_matrix.tsv")
   
   print(adata)

Loading mzTab Data
~~~~~~~~~~~~~~~~~~

You can download an example mzTab `file here <https://raw.githubusercontent.com/HUPO-PSI/mzTab/refs/heads/master/examples/1_0-Proteomics-Release/SILAC_SQ.mzTab>`_

.. code-block:: python

   # Load mzTab data
   adata = protdata.io.read_mztab("proteins.mzTab")
   
   print(adata)

Working with AnnData
~~~~~~~~~~~~~~~~~~~~

Once loaded, you can use all AnnData functionality:

.. code-block:: python

   # Basic operations
   adata.obs_names  # Sample names
   adata.var_names  # Protein names
   adata.X          # Main intensity matrix
   
   # Filtering
   adata = adata[adata.obs.index.str.contains("condition1"), :]  # Filter samples
   adata = adata[:, adata.var["Gene names"].notna()]             # Filter proteins
   
   # Save to h5ad format
   adata.write_h5ad("proteomics_data.h5ad")
   
   # Load from h5ad
   adata = anndata.read_h5ad("proteomics_data.h5ad")
