import React from 'react';
import { Group } from '@vx/group';
import { Tree } from '@vx/hierarchy';
import { LinearGradient } from '@vx/gradient';
import { hierarchy } from 'd3-hierarchy';

// import Links from './Links';
import Links from './LinksMove';

// import Nodes from './Nodes';
import Nodes from './NodesMove';
import {MenuItem, Select} from "@material-ui/core";

export default class extends React.Component {
  state = {
    layout: 'cartesian',
    orientation: 'horizontal',
    linkType: 'diagonal',
    stepPercent: 0.5
  };

  render() {
    const {
      data,
      width,
      height,
      events = false,
      margin = {
        top: 50,
        left: 100,
        right: 100,
        bottom: 50
      }
    } = this.props;
    const { layout, orientation, linkType, stepPercent } = this.state;

    if (width < 10) return null;

    const innerWidth = width - margin.left - margin.right;
    const innerHeight = height - margin.top - margin.bottom;

    let origin;
    let sizeWidth;
    let sizeHeight;

    if (layout === 'polar') {
      origin = {
        x: innerWidth / 2,
        y: innerHeight / 2
      }
      sizeWidth = 2 * Math.PI
      sizeHeight = Math.min(innerWidth, innerHeight) / 2
    } else {
      origin = { x: 0, y: 0 };
      if (orientation === 'vertical') {
        sizeWidth = innerWidth;
        sizeHeight = innerHeight;
      } else {
        sizeWidth = innerHeight;
        sizeHeight = innerWidth;
      }
    }

    const root = hierarchy(data, d => (d.isExpanded ? d.children : null))
    // root.each((node, i) => node.onClick = () => {
    //   console.log('clicked');
    // });

    return (
        <div>
          <div style={{"margin-top":"15px"}}>
            <label style={{"margin-left":"20px", "margin-right":"10px", "font-size": "1.2em", "font-family":"Helvetica", "font-weight":"400"}}>Layout:</label>
            <select style={{"margin-bottom":"10px", "border-radius":"5%", "border":"1px solid black", "font-size": "1.1em", "font-family":"Helvetica", "font-weight":"400"
              , "padding": "8px 16px", "background-color":"#FF6E90", "color":"white"}} onChange={e => this.setState({ layout: e.target.value })} value={layout}>
              <option value="cartesian">cartesian</option>
              <option value="polar">polar</option>
            </select>

            <label style={{"margin-left":"50px", "margin-right":"10px", "font-size": "1.2em", "font-family":"Helvetica", "font-weight":"400"}}>Orientation:</label>
            <select style={{"margin-bottom":"10px", "border-radius":"5%", "border":"1px solid black", "font-size": "1.1em", "font-family":"Helvetica", "font-weight":"400"
              , "padding": "8px 16px", "background-color":"#FF6E90", "color":"white"}}
                    onChange={e => this.setState({ orientation: e.target.value })} value={orientation} disabled={layout === 'polar'}>
              <option value="vertical">vertical</option>
              <option value="horizontal">horizontal</option>
            </select>
          </div>
          <svg width={width} height={height}>
            <LinearGradient id="lg" from="#fd9b93" to="#fe6e9e" />
            <rect width={width} height={height} rx={14} fill="#272b4d" />
            <Tree
                top={margin.top}
                left={margin.left}
                root={root}
                size={[
                  sizeWidth,
                  sizeHeight
                ]}
                separation={(a, b) => (a.parent == b.parent ? 1 : .5) / a.depth}
            >
              {({ data }) => (
                  <Group
                      top={origin.y}
                      left={origin.x}
                  >
                    <Links
                        links={data.links()}
                        linkType={linkType}
                        layout={layout}
                        orientation={orientation}
                        stepPercent={stepPercent}
                    />

                    <Nodes
                        nodes={data.descendants()}
                        layout={layout}
                        orientation={orientation}
                        onNodeClick={node => {
                          if (!node.data.isExpanded) {
                            node.data.x0 = node.x;
                            node.data.y0 = node.y;
                          }
                          node.data.isExpanded = !node.data.isExpanded;
                          this.forceUpdate()
                        }}
                    />

                  </Group>
              )}
            </Tree>
          </svg>
        </div>
    );
  }

}
