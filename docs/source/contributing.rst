Contributing
============

We welcome contributions to protdata! This guide will help you get started with contributing to the project.

Getting Started
---------------

1. Fork the repository on GitHub
2. Clone your fork locally:

.. code-block:: bash

   git clone https://github.com/yourusername/protdata.git
   cd protdata

3. Create a development environment:

.. code-block:: bash

   pip install -e ".[dev,docs]"

Development Workflow
--------------------

1. Make your changes
2. Run the tests:

.. code-block:: bash

   pytest tests/

Code Style
----------

- Follow PEP 8 style guidelines
- Use type hints where appropriate
- Add docstrings to all public functions
- Use numpy-style docstrings

Documentation
-------------

- Update docstrings for any changed functions
- Update tutorials if adding new features
- Build docs locally to check formatting:

.. code-block:: bash

   cd docs
   make html

Adding New Loaders
------------------

When adding support for a new proteomics tool:

1. Create a new loader function in `protdata/io/`
2. Follow the existing pattern:
   - Return an AnnData object
   - Use consistent variable naming (X, obs, var, uns, layers)
   - Add proper error handling
   - Include comprehensive docstrings

3. Add the loader to `protdata/io/__init__.py`
4. Write comprehensive tests
5. Add tutorial documentation
6. Update the API reference


Contact
-------

- GitHub Issues: `GitHub Issues <https://github.com/czbiohub-sf/protdata/issues>`_

Code of Conduct
---------------

This project follows the standard open-source code of conduct. Be respectful and inclusive in all interactions. 