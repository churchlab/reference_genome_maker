reference_genome_maker
======================

Python scripts to create reference genome files (e.g. Genbank) from various inputs (e.g. vcf)


### Requirements

We use virtualenv.

    $ virtualenv venv

    $ source venv/bin/activate .

    (venv) $ pip install -r requirements.txt

Do not modify `requirements.txt` by hand. After installing new requirements through pip, run:

    (venv) $ pip freeze > requirements.txt

### Testing

We use nose for testing. Running the following from the top-level directory should work:

    (venv) $ nosetests
