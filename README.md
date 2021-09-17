# PhyTReeV

Web app for reconstructing phylogenetic trees for different types of inputs using various algorithms.
This app is a part of my Master Thesis in Bioinformatics on Faculty of Mathematics, department of Informatics. 
You may find the thesis in Serbian.

Possible inputs would be:
- Distance Matrix
- Multiple Sequence Alignment
- Starting Tree For Small Parsimony

Depending on the input value we can then choose which algorithm we want to execute:
- UPGMA
- Neighbor-Joining
- Small Parsimony (Only with starting tree)
- Additive Phylogeny

For starting the BE you need to do the following steps (Linux, Mac):
- Position yourself inside the master-backend folder `cd master-backend`
- Download and install latest version of Python
- Install virtualenv `python3 -m pip install virtualenv`
- Creating and activating virtualenv 
`virtualenv .
source bin/activate`
- Installing Flask `python3 -m pip install Flask`
- Running the Flask server `flask run`

For starting the FE part of the app you need to do the following steps (Linux, Mac):
- Position yourself inside the master-backend folder `cd master-frontend`
- You'll need to install Node >= 14.0.0 and npm >= 5.6 https://nodejs.org/en/
- Then you should install Yarn with npm `npm install --global yarn`
- Run `yarn install` inside directory
- Run `yarn start`

Voilà, your app is up and running!
