function monteCarloPick(startNode, endNode){
    let siteseeds = [];
    /*
    seed start and end nodes // establish direction vector
    while randomizing course direction on random walk within radial range
    projection on course direction vector axis establishes weighting 
    for course correction.  
    */    

    ellipse(endNode.x,endNode.y,1,1);
    let rdir = 90*(2*PI/360);
    let walknodes = [startNode]
    var pick = createVector(startNode.x,startNode.y);
    var end = createVector(endNode.x, endNode.y);
    var distV = pick.sub(end);
    var distElapse = distV.mag();
    var t = 0;
    while (distElapse > 10 && t < 10000){
        let ppoint = walknodes[walknodes.length-1];
        pick = createVector(ppoint.x,ppoint.y);
        end = createVector(endNode.x,endNode.y);
        let npdir = rdir*random(.1,.7)*random([-1,1]);
        let distV = end.copy().sub(pick.copy()).rotate(npdir);
        distV.normalize().mult(random(8,16));
        let npick = pick.copy().add(distV);
        //console.log(npick);
        let npend = npick.copy().sub(end.copy());
        distElapse = npend.mag();
        walknodes.push({x:npick.x, y:npick.y});
        t+=1;
    }
    walknodes.push(endNode);
    //console.log(walknodes);
    return walknodes;
  }

  function pickBnodes(nodes){
    let bnodesset = [];
    let start = nodes[0];
    let end = nodes[nodes.length-1];
    let svec = createVector(start.x,start.y);
    let evec = createVector(end.x,end.y);
    let dvec = evec.copy().sub(svec.copy());
    dvec.rotate(20*random([-1,1])*2*PI/360).mult(.6);
    let picki = 0;
    while (picki < nodes.length-1){
      let bnodeg = {};
      let apick = random([1,2,3,4,5,6,7,8,9,10,11,12]);
      let npick = picki + apick;
      if (npick > nodes.length-1){
        break;
      }
      let pstart = nodes[npick];
      let psvec = createVector(pstart.x,pstart.y);
      let pevec = psvec.add(dvec);
      let pend = {x:pevec.x, y: pevec.y};
      bnodeg.start = pstart;
      bnodeg.end = pend;
      bnodesset.push(bnodeg);
      picki = npick;
    }
    return bnodesset;
  }

  function setup() { 
    createCanvas(750, 500);
    background(0);
    noStroke();
    fill(250);
    let startNode = {x:random()*750, y: random()*500};
    let endNode = {x: random()*750, y: random()*500};
    let branchDir = endNode.copy().sub(startNode.copy())
    let wnodes = monteCarloPick(startNode,endNode);

    for (var i = 0; i < wnodes.length-4; i++){
      fill(250);
      noStroke();
      ellipse(wnodes[i].x, wnodes[i].y,1,1);
      noFill();
      stroke(255);
      strokeWeight(.5);
      let p1 = wnodes[i];
      let p2 = wnodes[i+1];
      //line(p1.x,p1.y,p2.x,p2.y);
       let p3 = wnodes[i+2];
       let p4 = wnodes[i+3];
      curve(p1.x,p1.y,p2.x,p2.y,p3.x,p3.y,p4.x,p4.y);
    }
  }

  function draw() { 
    // background(0);
   }