function storePolyPts(polyobj, polyverts, pt){
    let ptstrng = pt.x.toString()+","+pt.y.toString();
    if (ptstrng in polyobj){
      return false;
    }
    else{
      polyobj[ptstrng] = pt;
      polyverts.push([pt.x,pt.y]);
      return true;
    }
  }
  
  function getPolyPts(polyobj){
    pts = [];
    for (let key in polyobj){
      pt = [polyobj[key].x,polyobj[key].y];
      pts.push(pt);
    }
    pts.push(pts[0]);
    return [pts];
  }

  //store vertex string keyed to pair vertex with height rank for polygon area
  function buildVertexObj(vertexObj, pt, polyrank){
    let ptstrng = pt[0].toString()+","+pt[1].toString();
    if (ptstrng in vertexObj){
      vertexObj[ptstrng]+= polyrank;
      
    }
    else{
      vertexObj[ptstrng] = polyrank;
    }
  }

  function setPolyRank(vertexObj, poly, polyrank){
    for (var i = 0; i < poly.length-1; i++){
       buildVertexObj(vertexObj, poly[i], polyrank);
    }
  }

  function reScalePoly(poly){
    let npoly = poly[0];
    let rpoly = [];
    for (var i = 0; i < npoly.length; i++){
      let pt = npoly[i];
      let ptvec = createVector(pt[0],pt[1]);
      let plen = ptvec.mag();
      ptvec.normalize();
      ptvec.mult(1.0*plen);
      rpoly.push([ptvec.x, ptvec.y]);
    }
    return [rpoly];
  }
  
  function buildBoundingBox(poly){
    let xvals = [];
    let yvals = [];
    let polypts = poly[0];
    for (var i = 0; i < polypts.length; i++){
      xvals.push(polypts[i][0]);
      yvals.push(polypts[i][1]);
    }
    let minx = min(xvals);
    let maxx = max(xvals);
    let miny = min(yvals);
    let maxy = max(yvals);
    bboxpts = [];
    for (var i = floor(minx); i <= maxx; i++){
      for (var j = floor(miny); j <= maxy; j++){
        bboxpts.push([i,j]);
      }
    }
    return bboxpts;
  }
  
  function getPointsInPoly(poly){
    //poly input is in turf js format
    let bpts = buildBoundingBox(poly);
    let points = turf.points(bpts);
    var searchWithin = turf.polygon(poly);
    var ptsWithin = turf.pointsWithinPolygon(points, searchWithin);
    return ptsWithin;
  }
  
  function pointOnVector(u, m, p1, p2){
    let pos = m.mult(u);
    let n1 = p1.add(pos);
    let t1 = p1.x >= n1.x && n1.x <= p2.x;
    let t2 = p1.y >= n1.y && n1.y <= p2.y;
    if (t1 && t2){
      return {n1: n1, check:true};
    }
    return {n1: n1, check:false};
  }

  function computeDenom(p1, p2, p3, p4){
    return ((p1.x - p2.x)*(p3.y-p4.y)-(p1.y-p2.y)*(p3.x-p4.x));
    //bezier line to line intersection
  }



  function computeNumer(p1,p2,p3,p4){
    return ((p1.x-p3.x)*(p3.y-p4.y)-(p1.y-p3.y)*(p3.x-p4.x));
  }

  function computeNumeru(p1,p2,p3,p4){
    return ((p2.x-p1.x)*(p1.y-p3.y)-(p2.y-p1.y)*(p1.x-p3.x));
  }

  function buildProfile(){
    let heightprofile = [0,.3,.35,.45,.55,.7,.75,.97,.98,.998,1];
    /*
    let profilelen = random()*20;
    for (var i = 0; i < profilelen; i++){
      heightprofile.push(random());
    }
    heightprofile.push(1.0);
    heightprofile.sort(function(a, b){return a-b});
    */
    return heightprofile;
  }

  function getProfileheight(heightprofile, t){
    let inc = 1/heightprofile.length;
    //floor(t) >= index *inc  means floor(t) /inc <= index
    let index = floor(t /inc)-1;
    let index2 = floor(t / inc);
    let hpmax = heightprofile[index2];
    let hpmin = heightprofile[index];
    return lerp(hpmin,hpmax, t);
  }
  
  function dScale2(centroid, poly, pt){
    let dists = [];
    var p3 = createVector(centroid[0], centroid[1]);
    var npt = createVector(pt[0],pt[1]);
    //let n = npt.sub(p3);
    //n.normalize();
    let ts = [];
    let maxdists = [];
    for (var i = 0; i < poly.length-1; i++){
      let p1 = createVector(poly[i][0], poly[i][1]);
      let p1p3 = p3.sub(p1);
      maxdists.push(p1p3.mag());

    }
    let maxdist = max(maxdists);
    //console.log(maxdist);
    p3 = createVector(centroid[0], centroid[1]);
    let nptp3 = npt.sub(p3);
    nptp3.normalize();
    nptp3.mult(1.5*maxdist);
    npt = p3.add(nptp3).copy();
    //console.log(npt);
    p3 = createVector(centroid[0], centroid[1]);
    
    //console.log(p3);
    for (var i = 0; i < poly.length-1; i++){
      let coord = poly[i];
      let coord2 = poly[i+1];
      
      let p1 = createVector(coord[0], coord[1]);
      let p2 = createVector(coord2[0], coord2[1]);
      let dres = computeDenom(p1,p2,p3,npt);
      if (abs(dres) == 0.000){
        console.log ('zero division error.');
      }
      let nres = computeNumer(p1,p2,p3,npt);
      let nres2 = computeNumeru(p1,p2,p3,npt);
      let t = nres/dres;
      let u = nres2/dres;
      let test = 0 <= t && t <= 1;
      let test2 = u >= 0 && u <= 1;
      ts.push(t);
      if (test && test2){
        let m = p2.sub(p1);
        let mt = m.mult(t);
        let npt2 = p1.add(mt);
        let npt2c = npt2.copy();
        let nptp3 = npt2.sub(p3);
        let lowaddhf = getVorVertexHeight(vertexObj, coord, coord2, [npt2c.x,npt2c.y]);
        return [nptp3.mag(),lowaddhf];
      }
    //   let m = p2.sub(p1);
    //  // m.normalize();
    //   // let p3 = createVector(centroid[0], centroid[1]);
    //   // let npt = createVector(pt[0],pt[1]);
    //   // let n = npt.sub(p3);
    //   // n.normalize();
    //   let marray = [[m.x, -n.x], [m.y, -n.y]];
    //   let mmatrix = math.matrix(marray);
    //   if (math.det(mmatrix) == 0){
    //     console.log('not invertible!');
    //     continue;
    //   }
    //   let nmp1 = p3.sub(p1);
    //   let nmp1arr = [[nmp1.x],[nmp1.y]];
    //   let nmp1matrix = math.matrix(nmp1arr);
    //   let invmmatrix = math.inv(mmatrix);
    //  // console.log (invmmatrix);
    //   //console.log (nmp1matrix);
    //   let sol = math.multiply(invmmatrix,nmp1matrix);
    //   let u = math.subset(sol, math.index(0, 0));

      //console.log(u);
      // let res = pointOnVector(u, m, p1, p2);
      // if (res.check){
      //   let n1 = res.n1;
      //   let p3n1 = p3.sub(n1);
      //   let p3n1d = p3n1.mag();
      //   dists.push(p3n1d);
      //   return p3n1d;
        
      // }
    }
    //console.log(ts);
    console.log('no intersection found!');
    return null;
  }
  
  function dScale(centroid, poly){
    let dists = [];
    for (var i = 0; i < poly.length-1; i++){
      let coord = poly[i];
      let coord2 = poly[i+1];
      
      let p1 = createVector(coord[0], coord[1]);
      let p2 = createVector(coord2[0], coord2[1]);
      let m = p2.sub(p1);
      m.normalize();
      
      let p3 = createVector(centroid[0], centroid[1]);
     
      let pma = p3.sub(p1);
      let pmadotm = pma.dot(m);
      let pmadotmmultm = m.mult(pmadotm);
      let d = pma.sub(pmadotmmultm).mag();
      //let tnot = (m.dot(pma))/(m.dot(m));
      //let tnotm = m.mult(tnot);
      //let pintrsct = p1.add(tnotm);
      //let p4 = pintrsct.sub(p3);
      //let p4d = p4.mag();
      dists.push(d);
    }
    return min(dists);
  }

  function getVaval(aval){
    let range1 = .01
    let range2 = .0415;
    let range3 = .0725;
    let range4 = .1415;
    let range5 = .3515;
    let range6 = .4615;
    let range7 = .5715;
    let range8 = .7815;
    let range9 = .915;
    if (aval >= range1 && aval <= range2){
      return lerp(.75, .85, aval);
    }    
    else if (aval > range2 && aval <= range3){
      return lerp(.85,.55, aval);
    }
    else if (aval > range3 && aval <= range4){
      return lerp(.55, .43, aval);
    }
    else if (aval > range4 && aval <= range5){
      return lerp(.43, .35, aval);
    }
    else if (aval > range5 && aval <= range6){
      return lerp(.35, .2715, aval);
    }
    else if (aval > range7 && aval <= range8){
      return lerp(.2715, .25, aval);
    }
    else if (aval > range8 && aval <= range9){
      return lerp(.25, .1, aval);
    }
    else if (aval > range9){
      return lerp(.1, 0, aval);
    }
    else if (aval < range1){
      return lerp(.41,.75, aval);
    }
  }

  function getVorVertexHeight(vertexObj, pt1, pt2, ipt){
    //ipt is point on line segment between pt1, and pt2

    let ptstrng = pt1[0].toString()+","+pt1[1].toString();
    let ptstrng2 = pt2[0].toString()+","+pt2[1].toString();
    //console.log(ptstrng);
    //console.log(ptstrng2);
    let aval = vertexObj[ptstrng];
    let aval2 = vertexObj[ptstrng2];
    //console.log(aval);
    //console.log(aval2);
    let pt1hf = getVaval(aval);
    let pt2hf = getVaval(aval2);
    let pt1vec = createVector(pt1[0],pt1[1]);
    let pt1vecc = pt1vec.copy();
    let pt2vec = createVector(pt2[0],pt2[1]);
    let iptvec = createVector(ipt[0],ipt[1]);
    let p12d = pt1vec.sub(pt2vec).mag();
    let p1id = pt1vecc.sub(iptvec).mag();
    let t = p1id/p12d;
    if (p12d == 0){
      console.log('zero division error getVorVertexHeight t param!');
    }
    //console.log(t);
    return lerp(pt1hf, pt2hf, t);
  }

  function getHeightFactor(){
    let clen = currenta-mina;
    let range1 = .025;
    //let rangemax = .4;
    let range2 = .05;
    let range3 = .1;
    let range4 = .2;
    let aval = clen/normlen;
    if (aval >= range1 && aval <= range2){
      let ns = range2-range1;
      let nsval = aval - range1;
      let t = nsval/ns;
      let hf = lerp(.25, 1, t);
      return hf;
    }
    else if ( aval > range2 && aval <= range3){
      let ns = range3-range2;
      let nsval = aval - range2;
      let t = nsval/ns;
      let hf = lerp(1, .5, t);
      return hf;      
    }

    else if (aval > range3 && aval <= range4){
      let ns = range4-range3;
      let nsval = aval - range3;
      let t = nsval/ns;
      let hf = lerp(.5, .05, t);
      return hf;      
    }
    else if (aval > range4){
      let ns = 1-range4;
      let nsval = aval - range4;
      let t = nsval/ns;
      let hf = lerp(.05, 0, t);
      return hf;      
    }
    else if (aval < range1){
      let ns = range1;
      let nsval = range1 - aval;
      let t = nsval/ns;
      let hf = lerp(0, .25, t);
      return hf;  
    }

    return clen/normlen;
  }

  function overlay(a,b){
    if (a < 0.5){
      return 2*a*b;
    }
    else{
      return (1 - 2*(1-a)*(1-b));
    }
  }
  
  function computeHeight2(centroid, interiorP,poly){
    //centroid of the polygon
    // interior points of the polygon (pixel array) turf js format
    let features = interiorP.features;
    let scaleD2 = dScale(centroid, poly);
    let scaleDset = [];
    let heightprofile = buildProfile();
    for (var i = 0; i < features.length; i++){
      let geometry = features[i].geometry;
      //console.log(geometry);
      let coord = geometry.coordinates;
      let coordstrng = coord[0].toString()+","+coord[1].toString();
      if (coordstrng in completedPoints){
        console.log('found repetition!');
        continue;
      }
      let p1 = createVector(centroid[0],centroid[1]);
      let p2 = createVector(coord[0],coord[1]);
      var dsreturn = dScale2(centroid, poly, coord);
      let scaleD = dsreturn[0];
      let lowaddhf = dsreturn[1];
      //console.log(lowaddhf);
      if (scaleD < 0){
        console.log(scaleD);
      }
      if (scaleD == null){
        console.log('hit');
        scaleD = scaleD2;
      }
      if (abs(scaleD) === 0.000){
        console.log('Zero Division ScaleD!');
      }
      
      let p3 = p2.sub(p1);
      let p3d = p3.mag();
      scaleDset.push(p3d/scaleD);
      let nh = getProfileheight(heightprofile,(p3d/scaleD));
      let hf = getHeightFactor();
      //console.log(lowaddhf);
      let nh1 = nh;
      var hf1 = lowaddhf*((1-min([1,nh1]))+1)/2;
      let hf2 = noise(coord[0]*.008,coord[1]*.008);
      let hf2b = noise(coord[0]*.1,coord[1]*.1);
      let hf2a = lerp(hf1,hf2,0.5);
      let hf3 = overlay(hf1,hf2);
      let hf3a = overlay(hf2,hf1);
      let hf3b = overlay(hf1,hf2b);
      let hf3c = overlay(hf2b,hf1);
      let hf4 = lerp(hf3,hf3a,0.5);
      let hf4a = lerp(hf3b,hf3c,.5);
      let hf5 = lerp(hf3,hf3b,.5);
      // if (hf1 >= .7){
      //   hf1 = 0;
      // }
      let hpt = hf1*255;

      //stroke(hpt);
      fill(hpt);
      noStroke();
      ellipse(coord[0],coord[1],1,1);
      completedPoints[coordstrng] = 1;
    }
    //console.log(scaleDset);
    stroke(255);
  }
  
  function computeHeight(centroid, interiorP,poly){
    //centroid of the polygon
    // interior points of the polygon (pixel array) turf js format
    let features = interiorP.features;
    let scaleD = dScale(centroid, poly);
    for (var i = 0; i < features.length; i++){
      let geometry = features[i].geometry;
      //console.log(geometry);
      let coord = geometry.coordinates;
      let p1 = createVector(centroid[0],centroid[1]);
      let p2 = createVector(coord[0],coord[1]);
      let p3 = p2.sub(p1);
      let p3d = p3.mag();
      let hpt = (1-min([1,p3d/scaleD]))*255;
      fill(hpt,10);
      noStroke();
      //strokeWeight(.0001);
      ellipse(coord[0],coord[1],2,2);
    }
    stroke(255);
  }

  var normlen;
  var mina;
  var currenta;
  var vertexObj;
  var completedPoints;
  
  function setup() { 
    createCanvas(750, 500);
    background(0);
    completedPoints = {};
    ///*
    let snum = 120000;
    let siteseedsnum = random(.5,1)*700;
    let siteseeds = [];
    var sitesobj = {};
    for (var i = 0; i < siteseedsnum; i++){
      let pick1y = random([-1,1]);
      var pick2y = 0;
      let pick1x = random([-1,1]);
      var pick2x = 0;
      if (pick1y == -1){
        pick2y = 1;
      }
      if (pick1x == -1){
        pick2x = 1;
      }
      siteseeds.push({x: (pick2x + pick1x*random(random()))*839, y:(pick2y+pick1y*random(1-random()))*639});
    }
    var sites = [];//{x:300,y:300}, {x:100,y:100}, {x:200,y:500}, {x:250,y:450}, {x:600,y:150}];
    for (var i = 0; i < snum; i++){
      let seedpick = random(siteseeds);
      let vlen = random(random())*30;
      let vrot = 2*PI*random(random());
      let v = createVector(0,1);
      v.rotate(vrot);
      v.mult(vlen);
      let a = createVector(seedpick.x,seedpick.y);
      a.add(v);
      let astrng = a.x.toString()+","+a.y.toString();
      if (astrng in sitesobj){
        continue;
      }
      sites.push({x:a.x,y:a.y});
      sitesobj[astrng] = 1;
    }
    console.log(sites);
    // xl, xr means x left, x right
    // yt, yb means y top, y bottom
    var bbox = {xl:-100, xr:900, yt:-100, yb:700};
    var voronoi = new Voronoi();
    // pass an object which exhibits xl, xr, yt, yb properties. The bounding
    // box will be used to connect unbound edges, and to close open cells
    background(0);
    result = voronoi.compute(sites, bbox);
    let edgs = result.edges;
    let cells = result.cells;
    var polys = [];
    let pareas = [];
    vertexObj = {};
    for (var i =0; i < cells.length; i++){
      let poly = {};
      let polyverts = [];
      let cell = cells[i];
      if (i == 0){
        console.log(cell);
      }
      let halfedges = cell.halfedges;
      for (var j = 0; j < halfedges.length; j++){
        let edge = halfedges[j].edge;
  
        let pt1 = halfedges[j].getStartpoint();
        let pt2 = halfedges[j].getEndpoint();
        //console.log(pt1);
        //console.log(pt2);
        storePolyPts(poly, polyverts,pt1);
        storePolyPts(poly, polyverts,pt2);
      }
      //polys.push(poly);
      //console.log(cell);
      ///*
      let poly1 = getPolyPts(poly);
      //console.log(poly1);
      var polygon = turf.polygon(poly1);
      var area = turf.area(polygon);
      setPolyRank(vertexObj, poly1[0], area);
      pareas.push(area);

    }
    mina = min(pareas);
    var maxa = max(pareas);
    normlen = maxa-mina;
    console.log('minimum area:')
    console.log(mina);
    console.log(maxa);
    for (let key in vertexObj){
      
      //console.log((vertexObj[key] - 3*mina)/(3*normlen));
      vertexObj[key] = (vertexObj[key] - 3*mina)/(3*normlen);
    }
    //console.log(vertexObj);
    polys = [];
    for (var i =0; i < cells.length; i++){
      let poly = {};
      let polyverts = [];
      let cell = cells[i];
      currenta = pareas[i];
      if (i == 0){
        console.log(cell);
      }
      let halfedges = cell.halfedges;
      for (var j = 0; j < halfedges.length; j++){
        let edge = halfedges[j].edge;
  
        let pt1 = halfedges[j].getStartpoint();
        let pt2 = halfedges[j].getEndpoint();
        //console.log(pt1);
        //console.log(pt2);
        storePolyPts(poly, polyverts,pt1);
        storePolyPts(poly, polyverts,pt2);
      }
      polys.push(poly);
      //console.log(cell);
      ///*
      let poly1 = getPolyPts(poly);
      let poly2 = reScalePoly(poly1);
      //console.log(poly2);
      var polygon = turf.polygon(poly1);
      var centroid = turf.centroid(polygon);
      let coord = centroid.geometry.coordinates;
      //if (i == 0){
        //console.log("points: ");
        //console.log(getPointsInPoly(poly1));
        ///*
        let interiorP = getPointsInPoly(poly1);
        //let scaleD = dScale(centroid, poly1[0]);
        computeHeight2(coord, interiorP, poly1[0]);
        
      //}
  
      
      //console.log(coord);
      //ellipse(coord[0], coord[1], 3, 3);
      
    }
    console.log(polys[0]);
    let poly = getPolyPts(polys[0]);
    console.log(poly);
    var polygon = turf.polygon(poly);
    var centroid = turf.centroid(polygon);
    console.log(centroid.geometry.coordinates);
    for (var i = 0; i < 750; i++){
      for (var j = 0; j < 500; j++){
        let ptstrng = i.toString()+","+j.toString();
        if (ptstrng in completedPoints){
          continue;
        }
        console.log('missing point');
        stroke(0);
        ellipse(i,j,.01,.01);
      }
    }
    // stroke(250);
    // for (var i = 0; i < edgs.length;i++){
    //   let edg = edgs[i];
    //   //line(edg.va.x,edg.va.y,edg.vb.x,edg.vb.y);
    //   //ellipse(edg.lSite.x, edg.lSite.y,2,2);
    // }
    // render, further analyze, etc.
    //*/
  } 
  
  function draw() { 
   // background(0);
  }