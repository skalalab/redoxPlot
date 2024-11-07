# Redox Project

## Visualizations

### Installation 
- Main dependencies: 
    - shiny, shinylive, matplotlib, pandas, plotly

### Deployment
- Use `shinylive export . docs  --full-shinylive` to deploy the shiny app
- Use `git lfs track "*.json"` to track large json files produced by shinylive so that they can be pushed to github. 
- `.github/workflows/static.yml` sets the directory for deployed/static assets and directs github to look for the actual lfs assets (json files), instead of its pointers.


# Update 
- now it is also deployed on [redoxPlot](https://allanware.shinyapps.io/redoxplot/)
    - shinylive uses pyodide which does not have the dash_bio package pre-included. I haven't figured out how to include it myself. 