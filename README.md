# Redox Project

## Visualizations

### Installation 
- See *requirements.txt* for dependencies. Main ones: 
    - shiny, shinylive, matplotlib, pandas

### Deployment
- Use `shinylive export . docs  --full-shinylive` to deploy the shiny app
- Use `git lfs track "*.json"` to track large json files produced by shinylive so that they can be pushed to github. 
- `.github/workflows/static.yml` sets the directory for deployed/static assets and directs github to look for the actual lfs assets (json files), instead of its pointers. 