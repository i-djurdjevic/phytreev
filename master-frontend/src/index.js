import React from 'react';
import {render} from 'react-dom';
import Tree from './tree/Tree';
import data from './data';
import AppBar from '@material-ui/core/AppBar';
import CssBaseline from '@material-ui/core/CssBaseline';
import Drawer from '@material-ui/core/Drawer';
import TextField from '@material-ui/core/TextField';
import Hidden from '@material-ui/core/Hidden';
import IconButton from '@material-ui/core/IconButton';
import List from '@material-ui/core/List';
import ListItem from '@material-ui/core/ListItem';
import ListItemIcon from '@material-ui/core/ListItemIcon';
import ListItemText from '@material-ui/core/ListItemText';
import LoopIcon from '@material-ui/icons/Loop';
import MenuIcon from '@material-ui/icons/Menu';
import Toolbar from '@material-ui/core/Toolbar';
import Typography from '@material-ui/core/Typography';
import {makeStyles, useTheme} from '@material-ui/core/styles';
import {applyAlgorithm} from "./algorithmsService";
import {tree_newick_parser} from "./utils/treeNewickParser";
import {
  parseAdditivePhylogenySteps,
  parseArrayTable,
  parseInitialMatrix,
  parseInitialParsimonySteps,
  parseNjMatrix,
  parseNjSteps,
  parseSmallParsimonySteps,
  parseTable,
  parseUpgmaSteps
} from "./utils/stepsParser";
import Radio from '@material-ui/core/Radio';
import RadioGroup from '@material-ui/core/RadioGroup';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import {Button} from "@material-ui/core";
// import { Matrix } from './matrix';

const drawerWidth = 400;

const useStyles = makeStyles((theme) => ({
  main: {
    marginBottom: '100px'
  },
  root: {
    display: 'flex',
  },
  drawer: {
    [theme.breakpoints.up('sm')]: {
      width: drawerWidth,
      flexShrink: 0,
    },
  },
  appBar: {
    [theme.breakpoints.up('sm')]: {
      width: `calc(100% - ${drawerWidth}px)`,
      marginLeft: drawerWidth,
    },
  },
  menuButton: {
    marginRight: theme.spacing(2),
    [theme.breakpoints.up('sm')]: {
      display: 'none',
    },
  },
  toolbar: theme.mixins.toolbar,
  drawerPaper: {
    width: drawerWidth,
  },
  content: {
    flexGrow: 1,
  },
  textField: {
    position: 'relative',
    width: '350px',
    left: '5%',
    top: '5%',
    alignContent: 'center'
  },
  board: {
    float: 'left'
  },
  radioButtons: {
    position: 'relative',
    left: '5%',
  },
  algorithmList: {
    position: 'relative',
    top: '5%',
  },
  stepsTextField: {
    position: 'relative',
    top: '5%',
    left: '5%',
    width: '350px',
    alignContent: 'center',
    marginBottom: '10px'
  },
  stepsTableField: {
    position: 'relative',
    top: '5%',
    marginBottom: '20px',
    width: '350px',
    margin: 'auto'
  },
  fileUpload: {
    position: 'relative',
    left: '5%',
    top: '2%'
  },
  input: {
    display: 'none'
  },
  stepsButton: {
    margin: '10px'
  },
  table: {
    minWidth: 650,
  },
}));

let tree = tree_newick_parser("(((2,3)5,1)4)0;");
const initialMatrix =
    "0 30 34 21\n"
    + "30 0 28 39\n"
    + "34 28 0 43\n"
    + "21 39 43 0\n";

const App = (props) => {
  const {window} = props;
  const classes = useStyles();
  const theme = useTheme();
  const [mobileOpen, setMobileOpen] = React.useState(false);
  const [matrix, setMatrix] = React.useState(initialMatrix);
  const [treeData, setTreeData] = React.useState(tree);
  const [radio, setRadioValue] = React.useState('distance');
  const [textFieldLabel, setTextFieldLabel] = React.useState('Distance Matrix');
  const [stepsTextFieldLabel, setStepsTextFieldLabel] = React.useState(
      'Choose Algorithm');
  const [stepsValue, setStepsValue] = React.useState('');
  const [numberOfSteps, setNumberOfSteps] = React.useState([]);
  const isMountedRef = React.useRef(true)
  const [selectedFile, setSelectedFile] = React.useState(null);
  const [isFilePicked, setIsFilePicked] = React.useState(false);
  const [uploadText, setUploadText] = React.useState("Choose file");
  const [result, setResult] = React.useState();
  const [stepTable0, setStepTable0] = React.useState();
  const [stepTable1, setStepTable1] = React.useState();
  const [stepTable2, setStepTable2] = React.useState();
  const [algorithm, setAlgorithm] = React.useState("");

  const handleDrawerToggle = () => {
    setMobileOpen(!mobileOpen);
  };

  React.useEffect(() => () => {
    isMountedRef.current = false
  }, [])

  const resetAllValues = () => {
    setTreeData(data);
    setStepTable0(null);
    setStepTable2(null);
    setStepTable1(null);
    setStepsValue("");
    setStepsTextFieldLabel("Choose Algorithm");
    setNumberOfSteps([]);
  }

  function onClickAlgorithm(inputMatrix, alg) {
    setStepTable1(null);
    if (alg === "Neighbor Joining") {
      applyAlgorithm("nj", inputMatrix, radio)
      .then((result) => {
            if (isMountedRef.current) {
              setAlgorithm("nj");
              setResult(result);
              setNumberOfSteps(result.steps);
              setStepsTextFieldLabel('Neighbor Joining');
              setStepsValue(parseNjSteps(result, 0));
              setStepTable2(parseNjMatrix(result, result.dStarMatrices, 0));
              setStepTable1(parseTable(result, 0));
              if (textFieldLabel !== "Multiple Sequence Alignment")
                setStepTable0(parseArrayTable(parseInitialMatrix(inputMatrix, result),
                  "Matrix:"));
              setTreeData(tree_newick_parser(result.tree));
            } else {
              console.log("not mounted");
            }
          }
      )
      .catch(() => {
        resetAllValues();
      });
    } else if (alg === "UPGMA") {
      applyAlgorithm("upgma", inputMatrix, radio)
      .then((result) => {
            if (isMountedRef.current) {
              setAlgorithm("upgma");
              setResult(result);
              setStepsTextFieldLabel('UPGMA');
              setStepsValue(parseUpgmaSteps(result, 0));
              if (textFieldLabel !== "Multiple Sequence Alignment")
                setStepTable0(parseArrayTable(parseInitialMatrix(inputMatrix, result),
                  "Matrix:"));
              setStepTable1(parseTable(result, 0));
              setStepTable2(null);
              setNumberOfSteps(result.steps);
              setTreeData(tree_newick_parser(result.matrix));
            } else {
              console.log("not mounted");
            }
          }
      )
      .catch(() => {
        resetAllValues();
      });
    } else if (alg === "Small Parsimony") {
      applyAlgorithm("small_parsimony", inputMatrix, radio)
      .then((result) => {
            if (isMountedRef.current) {
              setStepTable0(null);
              setStepTable2(null);
              setStepTable1(null);
              setResult(result);
              setTextFieldLabel("Small Parsimony");
              setAlgorithm("parsimony");
              setStepTable2(null);
              setTreeData(tree_newick_parser(result.tree));
              let parsedSteps = parseInitialParsimonySteps(result.steps);
              setNumberOfSteps(Object.keys(parsedSteps));
              setStepsValue(parseSmallParsimonySteps(result, parsedSteps,
                  parseInitialParsimonySteps(result.steps_calculations), 0));
            } else {
              console.log("not mounted");
            }
          }
      )
      .catch(() => {
        resetAllValues();
      });
    } else if (alg === "Additive Phylogeny") {
      applyAlgorithm("additive_phylogeny", inputMatrix, radio)
      .then((result) => {
            if (isMountedRef.current) {
              setResult(result);
              setStepsTextFieldLabel("Additive Phylogeny");
              setAlgorithm("additive")
              setTreeData(tree_newick_parser(result.tree));
              setStepsValue(parseAdditivePhylogenySteps(result, 0));
              if (result.limbLengthCalculations.length) {
                setNumberOfSteps(result.limbLengthCalculations);
              } else {
                setNumberOfSteps([]);
              }
              if (textFieldLabel !== "Multiple Sequence Alignment")
                setStepTable0(parseArrayTable(parseInitialMatrix(inputMatrix, result),
                  "Matrix:"));
              setStepTable1(parseArrayTable(result.dBalds[0], "D Bald:"));
              setStepTable2(parseArrayTable(result.dTrim[0], "D Trim:"));
            } else {
              console.log("not mounted");
            }
          }
      )
      .catch(() => {
        resetAllValues();
      });
    }
  }

  const handleTextFieldChange = (event) => {
    setMatrix(event.target.value);
    setUploadText("Choose file")
  };

  const handleStepButton = event => {
    let step = event.currentTarget.value;
    if (algorithm && algorithm === "upgma") {
      setStepsValue(parseUpgmaSteps(result, step));
      setStepTable1(parseTable(result, step));
    } else if (algorithm && algorithm === "nj") {
      setStepsValue(parseNjSteps(result, step));
      setStepTable2(parseNjMatrix(result, result.dStarMatrices, step));
      setStepTable1(parseTable(result, step));
    } else if (algorithm && algorithm === "additive") {
      setStepsValue(parseAdditivePhylogenySteps(result, step));
      setStepTable1(parseArrayTable(result.dBalds[step], "D Bald:"));
      setStepTable2(parseArrayTable(result.dTrim[step], "D Trim:"))
    } else if (algorithm && algorithm === "parsimony") {
      setStepsValue(parseSmallParsimonySteps(result,
          parseInitialParsimonySteps(result.steps),
          parseInitialParsimonySteps(result.steps_calculations), step));
    }
  }

  const handleRadioChange = event => {
    setRadioValue(event.target.value);
    if (event.target.value === "distance") {
      setTextFieldLabel("Distance Matrix");
    } else if (event.target.value === "alignment") {
      setTextFieldLabel("Multiple Sequence Alignment");
    } else {
      setTextFieldLabel("Starting tree for small parsimony");
    }
  };

  const changeHandler = event => {
    setSelectedFile(event.target.files[0]);
    setUploadText(event.target.files[0].name);
    setIsFilePicked(true);
  };

  const handleSubmission = () => {
    const formData = new FormData();
    formData.append('File', selectedFile);
    formData.append("type", radio);
    fetch('http://localhost:5000/upload_file',
        {
          method: 'POST',
          body: formData,
        }
    ).then((response) => response.json())
    .then((result) => {
      console.log('Success:', result);
      setMatrix(result.matrix);
      setUploadText("Data uploaded!")
    })
    .catch((error) => {
      console.error('Error:', error);
    });
  };

  const drawer = (
      <div className={classes.main}>
        <div className={classes.toolbar}/>
        <RadioGroup
            className={classes.radioButtons}
            aria-label="algorithmMethod" name="algorithmMethod"
            value={radio} onChange={handleRadioChange}>
          <FormControlLabel value="parsimony" control={<Radio/>}
                            label="Starting tree for small parsimony"/>
          <FormControlLabel value="alignment" control={<Radio/>}
                            label="Multiple sequence alignment"/>
          <FormControlLabel value="distance" control={<Radio/>}
                            label="Distance matrix"/>
        </RadioGroup>
        <div className={classes.fileUpload}>
          <Button variant="contained"
                  component="label">
            {uploadText}
            <input
                accept="csv*"
                onChange={changeHandler}
                hidden
                id="contained-button-file"
                multiple
                type="file"
            />
          </Button>
          <label style={{"margin-left": "20px"}} htmlFor="contained-button-file">
            <Button variant="contained" color="primary"
                    onClick={handleSubmission}>
              Upload
            </Button>
          </label>
        </div>
        <TextField
            fullWidth={true}
            id="matrix-textfield"
            className={classes.textField}
            label={textFieldLabel}
            multiline={true}
            margin="dense"
            rows={10}
            value={matrix}
            defaultValue={matrix}
            variant="outlined"
            type="string"
            onChange={handleTextFieldChange}
        />
        <List className={classes.algorithmList}>
          {['Small Parsimony', 'UPGMA',
            'Neighbor Joining', 'Additive Phylogeny'].map(
              (text, index) => (
                  <ListItem button key={text} onClick={() => {
                    onClickAlgorithm(matrix, text);
                  }}>
                    <ListItemIcon><LoopIcon/></ListItemIcon>
                    <ListItemText primary={text}/>
                  </ListItem>
              ))}
        </List>
        <TextField
            id="steps-textfield"
            className={classes.stepsTextField}
            label={stepsTextFieldLabel}
            multiline
            margin="dense"
            rows={10}
            variant="outlined"
            type="string"
            disabled={true}
            value={stepsValue}
        />
        <div className={classes.stepsTableField}>{stepTable0}</div>
        <div className={classes.stepsTableField}>{stepTable1}</div>
        <div className={classes.stepsTableField}>{stepTable2}</div>
      </div>
  );

  const container = window !== undefined ? () => window().document.body
      : undefined;

  return (
      <div className={classes.root}>
        <CssBaseline/>
        <AppBar position="fixed" className={classes.appBar}>
          <Toolbar>
            <IconButton
                color="inherit"
                aria-label="open drawer"
                edge="start"
                onClick={handleDrawerToggle}
                className={classes.menuButton}
            >
              <MenuIcon/>
            </IconButton>
            <Typography variant="h6" noWrap>
              Tree reconstruction and visualisation
            </Typography>
          </Toolbar>
        </AppBar>
        <nav className={classes.drawer} aria-label="mailbox folders">
          <Hidden smUp implementation="css">
            <Drawer
                container={container}
                variant="temporary"
                anchor={theme.direction === 'rtl' ? 'right' : 'left'}
                open={mobileOpen}
                onClose={handleDrawerToggle}
                classes={{
                  paper: classes.drawerPaper,
                }}
                ModalProps={{
                  keepMounted: true,
                }}
            >
              {drawer}
            </Drawer>
          </Hidden>
          <Hidden xsDown implementation="css">
            <Drawer
                classes={{
                  paper: classes.drawerPaper,
                }}
                variant="permanent"
                open
            >
              {drawer}
            </Drawer>
          </Hidden>
        </nav>
        <main className={classes.content}>
          <div className={classes.toolbar}/>
          <div className={classes.board}>
            <Tree data={treeData} width={document.documentElement.clientWidth
            - drawerWidth}
                  height={document.documentElement.clientHeight - 200}/>
          </div>
          {numberOfSteps.map((key, buttonValue) => {
            return (
                <Button className={classes.stepsButton} variant="contained"
                        color="secondary" onClick={e => handleStepButton(e)}
                        value={buttonValue}>{buttonValue}</Button>);
          })}
        </main>
      </div>
  );
};

render(<App/>, document.getElementById('root'));
