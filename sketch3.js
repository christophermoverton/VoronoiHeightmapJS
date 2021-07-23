function monteCarloPick(startNode, endNode){
    let siteseeds = [];
    /*
    seed start and end nodes // establish direction vector
    while randomizing course direction on random walk within radial range
    projection on course direction vector axis establishes weighting 
    for course correction.  
    */    
   function getCorrectionvector(npick){
     let nearpoints = kdt.nearest({x: npick.x, y: npick.y}, 50);
     //console.log(nearpoints);
     var corvec = createVector(0,0);
     var div = 0;
     for (var i = 0; i < nearpoints.length; i++){
       let npt = nearpoints[i];
       let nptvec = createVector(npt.x,npt.y);
       let nptnpick = npick.copy().sub(nptvec);
       if (nptnpick.mag() < 30){
         corvec.add(nptnpick);
         div += 1;
       }
     }
     if (div > 0){
       corvec.normalize().mult(-31);
       return corvec.copy();
     }
     return createVector(0,0);
   }

    ellipse(endNode.x,endNode.y,1,1);
    let rdir = 90*(2*PI/360);
    let walknodes = [startNode]
    var pick = createVector(startNode.x,startNode.y);
    var end = createVector(endNode.x, endNode.y);
    var distV = pick.sub(end);
    var distElapse = distV.mag();
    var t = 0;
    while (distElapse > 10 && t < 1000){
        let ppoint = walknodes[walknodes.length-1];
        pick = createVector(ppoint.x,ppoint.y);
        end = createVector(endNode.x,endNode.y);
        let npdir = rdir*random(.1,.7)*random([-1,1]);
        let distV = end.copy().sub(pick.copy()).rotate(npdir);
        distV.normalize().mult(random(2,5));
        var npick = pick.copy().add(distV);

        let corvec = getCorrectionvector(npick);
        npick.add(corvec);
        //console.log(npick);
        let npend = npick.copy().sub(end.copy());
        distElapse = npend.mag();
        walknodes.push({x:npick.x, y:npick.y});
        kdt.insert({x:npick.x, y:npick.y});
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

    let picki = 0;
    while (picki < nodes.length-1){
      let bnodeg = {};
      let apick = random([4,5,6,7,8,9,10]);
      let npick = picki + apick;
      if (npick > nodes.length-1){
        break;
      }
      let pstart = nodes[npick];
      let pnstart = nodes[npick-1];
      let psvec = createVector(pstart.x,pstart.y);
      let pnsvec = createVector(pnstart.x, pnstart.y);
      let dvec = svec.copy().sub(evec.copy());
      let dvec2len = svec.copy().sub(evec).mag();
      dvec.normalize();
      dvec.rotate(random(55,60)*random([-1,1])*2*PI/360).mult(random(.4,.1)*dvec2len);
      let pevec = psvec.add(dvec);
      let pend = {x:pevec.x, y: pevec.y};
      bnodeg.start = pstart;
      bnodeg.end = pend;
      bnodesset.push(bnodeg);
      picki = npick;
    }
    return bnodesset;
  }

  var kdt;

  function setup() { 
    createCanvas(750, 500);
    background(0);
    noStroke();
    fill(250);
    let startNode = {x:random(0,.1)*750, y: random()*500};
    let endNode = {x: random(.9,1)*750, y: random()*500};
    var distance = function(a, b){
      return Math.pow(a.x - b.x, 2) +  Math.pow(a.y - b.y, 2);
    }
    kdt = new kdTree([startNode],distance,["x","y"]);
    //let branchDir = endNode.copy().sub(startNode.copy())
    let wnodes = monteCarloPick(startNode,endNode);
    let bnodes = pickBnodes(wnodes);
    console.log(bnodes);
    for (var i = 0; i < wnodes.length-3; i++){
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

    for (var i = 0; i < bnodes.length; i++){
      let nstartNode = bnodes[i].start;
      let nendNode = bnodes[i].end;
      let newnodes = monteCarloPick(nstartNode,nendNode);
      console.log(newnodes);
      ///*
      for (var j = 0; j < newnodes.length-4; j++){
        fill(250);
        noStroke();
        ellipse(newnodes[j].x, newnodes[j].y,1,1);
        noFill();
        stroke(255);
        strokeWeight(.5);
        let p1 = newnodes[j];
        let p2 = newnodes[j+1];
        //line(p1.x,p1.y,p2.x,p2.y);
         let p3 = newnodes[j+2];
         let p4 = newnodes[j+3];
        curve(p1.x,p1.y,p2.x,p2.y,p3.x,p3.y,p4.x,p4.y);     
      }
      //*/
    }
  }

  function draw() { 
    // background(0);
   }