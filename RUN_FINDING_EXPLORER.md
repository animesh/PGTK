# Run the complete finding explorer

After a successful pipeline run:

```bash
cd results/igv/findings/finding_explorer
export PGTK_IGV_REPORTS_IMAGE=/cluster/projects/nn9036k/scrbkup/pgtk/singularity_cache/quay.io-biocontainers-igv-reports-1.16.0--pyh7e72e81_0.img
./serve_explorer.sh "$PWD" 8765
```

From the workstation, create an SSH tunnel:

```bash
ssh -L 8765:127.0.0.1:8765 ash022@login.saga.sigma2.no
```

Open `http://127.0.0.1:8765`. Selecting IGV generates and caches one standalone offline report for that finding. No finding is omitted from search or navigation.
