var nu_min  ;
var nu_max  ;
var nu_opt  ;

var numData;
var parameters;
var nu_i;
var N;
var N1;

var Pi = Math.PI;
var tau;
var n;
var h;
var kappa0 = [];  // curvature of the midline at time t = 0
var kappa = [];  // curvature of the midline at time t = 0

var plotClicked;

var init_parameters = function (nu_i) {
 
     numArray = numData[nu_i - 1];
    
  N = Math.round((numArray.length -3)/3);
   N1 = N - 1;
  h = 1 / N;

  // creating the init_parameters
};
var avg = function (arr) {
   return arr.reduce((a, b) => a + b, 0) / arr.length;
};

var hx0 = [];
var hy0 = [];
var hz0 = [];

var rx0 = [];
var ry0 = [];
var rz0 = [];



var cen = [0, 0, 0];

var data_t0 = function (nu_i) {
  h = 1/N;

  let tau = (1.29 + 0.015877763328999*(nu_i-1))*2*Pi;
  let tx0 = [];
  let ty0 = [];
  let tz0 = []; 

  let hx2p = [];
  let hy2p = [];
  let hz2p = [];

  kappa0 = Array(N + 1).fill(0);

    numArray = numData[nu_i - 1];
  // data read. now calculating quantities

 
  i = 0;
  hx0[i] = numArray[i * 3];
  hy0[i] = numArray[i * 3 + 1];
  hz0[i] = numArray[i * 3 + 2];

  for (var i = 1; i < N + 1; i++) {
    hx0[i] = numArray[i * 3];
    hy0[i] = numArray[i * 3 + 1];
    hz0[i] = numArray[i * 3 + 2];

    tx0[i] = (hy0[i - 1] * hz0[i] - hy0[i] * hz0[i - 1])/tau/h;
    ty0[i] = (hz0[i - 1] * hx0[i] - hz0[i] * hx0[i - 1])/tau/h;
    tz0[i] = (hx0[i - 1] * hy0[i] - hx0[i] * hy0[i - 1])/tau/h;

    

  }
  // if (isNaN(hx0  ) || isNaN(hy0  ) || isNaN(hz0   )) {
  //   console.log(`Found NaN at index ${i}`);
  // }

  tx0[0] = tx0[N];
  ty0[0] = ty0[N];
  tz0[0] = tz0[N];

  //let kmax = Math.max(...kappa0.map(Math.abs));;
 
// calculating curvature 
// curvature
hx2p[0] = (hx0[1] - hx0[N - 1] - 2 * hx0[0]) / h ** 2;
hy2p[0] = (hy0[1] - hy0[N - 1] - 2 * hy0[0]) / h ** 2;
hz2p[0] = (hz0[1] - hz0[N - 1] - 2 * hz0[0]) / h ** 2;

hx2p[N] = -hx2p[0];
hy2p[N] = -hy2p[0];;
hz2p[N] = -hz2p[0];;


 kappa0[0] = tx0[0] * hx2p[0] + ty0[0] * hy2p[0] + tz0[0] * hz2p[0];

for (var i = 1; i < N; i++) {
  hx2p[i] = (hx0[i + 1] + hx0[i - 1] - 2 * hx0[i]) / h ** 2;
  hy2p[i] = (hy0[i + 1] + hy0[i - 1] - 2 * hy0[i]) / h ** 2;
  hz2p[i] = (hz0[i + 1] + hz0[i - 1] - 2 * hz0[i]) / h ** 2;

  kappa0[i] = tx0[i] * hx2p[i] + ty0[i] * hy2p[i] + tz0[i] * hz2p[i];
}

 kappa0[N] = tx0[N] * hx2p[N] + ty0[N] * hy2p[N] + tz0[N] * hz2p[N];

//  // normalizing curvature 

   let kmax = Math.max(...kappa0.map(Math.abs));;
  kappa0 = kappa0.map((num) => num / tau);

 
///////////////////////


 
  // Midline
  rx0[0] = 0.0;
  ry0[0] = 0.0;
  rz0[0] = 0.0;
  h = 1;
  for (var i = 0; i < N; i++) {
    rx0[i + 1] = rx0[i] + h * tx0[i + 1];
    ry0[i + 1] = ry0[i] + h * ty0[i + 1];
    rz0[i + 1] = rz0[i] + h * tz0[i + 1];
  }

  cen[0] = avg(rx0);
  cen[1] = avg(ry0);
  cen[2] = avg(rz0);

  for (var i = 0; i < N + 1; i++) {
    rx0[i] = rx0[i] - cen[0];
    ry0[i] = ry0[i] - cen[1];
    rz0[i] = rz0[i] - cen[2];
  }
  // nomrlalizing the midline
  let dl = 0;
  for (let i = 0; i < N; i++) {
    let dx = rx0[i + 1] - rx0[i];
    let dy = ry0[i + 1] - ry0[i];
    let dz = rz0[i + 1] - rz0[i];
    dl += Math.sqrt(dx * dx + dy * dy + dz * dz);
  }
  dl = dl / 6;
  for (let i = 0; i <= N; i++) {
    rx0[i] /= dl;
    ry0[i] /= dl;
    rz0[i] /= dl;
  }

   
};

//  midline and binormal loaded from the saved file 

//the function below generates animation data 
var animationData = function(j) {
  let ind = j % (2*N); // ensuring that index is not out of bound
  let ind2 = j % (N); // ensuring that index is not out of bound
   if (ind == 0) {
     ind = 1;
  }
   let or = -1;

  if(ind<N+1){
 
  let ar1 = Array.from(hx0.slice(ind - 1, N));
  let ar2 = Array.from(hx0.slice(0, ind).map((x) => or * x));
  
  var hx = ar1.concat(ar2);

  ar1 = Array.from(hy0.slice(ind - 1, N ));
  ar2 = Array.from(hy0.slice(0, ind).map((x) => or * x));
  var hy = ar1.concat(ar2);

  ar1 = Array.from(hz0.slice(ind - 1, N));
  ar2 = Array.from(hz0.slice(0, ind).map((x) => or * x));
  var hz = ar1.concat(ar2);

  ar1 = Array.from(kappa0.slice(ind - 1, N ));
  ar2 = Array.from(kappa0.slice(0, ind).map((x) => or * x));
  var kappa = ar1.concat(ar2);

  // updated binormal 
   ar1 = Array.from(rx0.slice(ind - 1, N));
   ar2 = Array.from(rx0.slice(0, ind).map((x) => 1 * x));
  
   var rx = ar1.concat(ar2);



   ar1 = Array.from(ry0.slice(ind - 1, N));
   ar2 = Array.from(ry0.slice(0, ind).map((x) => 1 * x));
  
   var ry = ar1.concat(ar2);

   ar1 = Array.from(rz0.slice(ind - 1, N));
   ar2 = Array.from(rz0.slice(0, ind).map((x) => 1 * x));
  
   var rz = ar1.concat(ar2);
   
    
  }
  else {
  let ar1 = Array.from(hx0.slice(ind2-1, N).map((x) => or * x));
  let ar2 = Array.from(hx0.slice(0, ind2));

  var hx = ar1.concat(ar2);

   ar1 = Array.from(hy0.slice(ind2-1, N).map((x) => or * x));
   ar2 = Array.from(hy0.slice(0, ind2));

  var hy = ar1.concat(ar2);

   ar1 = Array.from(hz0.slice(ind2-1, N).map((x) => or * x));
   ar2 = Array.from(hz0.slice(0, ind2));

  var hz = ar1.concat(ar2);

  ar1 = Array.from(kappa0.slice(ind2-1, N).map((x) => or * x));
   ar2 = Array.from(kappa0.slice(0, ind2));

  var kappa = ar1.concat(ar2);

  // updated binormal 
   ar1 = Array.from(rx0.slice(ind2-1, N) );
   ar2 = Array.from(rx0.slice(0, ind2));

  var rx = ar1.concat(ar2);

   ar1 = Array.from(ry0.slice(ind2-1, N) );
   ar2 = Array.from(ry0.slice(0, ind2));

  var ry = ar1.concat(ar2);

   ar1 = Array.from(rz0.slice(ind2-1, N) );
   ar2 = Array.from(rz0.slice(0, ind2));

  var rz = ar1.concat(ar2);
 
  };


   // ensuring that centroid is at the origin 
   cen[0] = avg(rx);
  cen[1] = avg(ry);
  cen[2] = avg(rz);


  for (var i = 0; i < N + 1; i++) {
    rx[i] = rx[i] - cen[0];
    ry[i] = ry[i] - cen[1];
    rz[i] = rz[i] - cen[2];
  }
   /// rotating appropriately 
   
   // first rotating the hand such that index 0 is fixed
   
  let th =  -Math.acos(rx[0] / Math.sqrt(rx[0]*rx[0] + ry[0]*ry[0]));
   // angle hetween point at s==0 and x axis
   if (ry[0] < 0) {
    th = 2 * Math.PI - th;
  }
  
  if (ry[0] == 0 && rx[0] > 0) {
    th = 0;
  } else if (ry[0] == 0 && rx[0] < 0) {
    th = Math.PI;
  }
  //th = 0;
  norm1 =Math.sqrt(rx[0]**2 + ry[0] **2)
  //let costh = ry[0]/norm1;
  //let sinth = rx[0]/norm1;
  let costh = -rx[0]/norm1;
  let sinth = ry[0]/norm1;
  let tempx = 0;
  let tempy = 0;
   
    for (let i = 0; i <  N+1; i++) {
     
      
     tempx  = costh *hx[i]  - sinth *hy[i];
     tempy  = sinth *hx[i]  + costh *hy[i];
     hx[i] = tempx;
     hy[i] = tempy;

     tempx  = costh *rx[i]  - sinth *ry[i];
     tempy  = sinth *rx[i]  + costh *ry[i];
     rx[i] = tempx;
     ry[i] = tempy;
 

  };
  


 
  //console.log(rx ); // returns true if any element is NaN
  //console.log(rx.some(Number.isNaN)); // returns true if any element is NaN
  
   

  // cen[0] = avg(rx);
  // cen[1] = avg(ry);
  // cen[2] = avg(rz);
  
  // console.log(rx[240])
  //   console.log(rx.some(Number.isNaN));
  // for (var i = 0; i < N + 1; i++) {
  //   rx[i] = rx[i] - cen[0];
  //   ry[i] = ry[i] - cen[1];
  //   rz[i] = rz[i] - cen[2];
  // }

  
 
   
   
  return [hx,hy,hz,rx,ry,rz,kappa];
  // rotating appropriately
}
// matrix multiplication

 
 
var hx = [];
var hy = [];
var hz = [];

var rx = [];
var ry = [];
var rz = [];
var kappa = [];
var fourierExpansion = function (n, t, hl) {

   l = (hl * 1.7) / 2;
   var v = [];
  
   console.log(l);
  const result = animationData(t);
   hx = result[0];
   hy = result[1];
   hz = result[2];
  
   rx = result[3];
   ry = result[4];
   rz = result[5];

   for(var i=0; i<N+1;i++) {
    //rz[i]=rz[i]-.35;
    //rz[i]=rz[i]-.35;
   };

   kappa = result[6];
  
   cen[0] = avg(rx);
   cen[1] = avg(ry);
   cen[2] = avg(rz);
 
  var mid = [0, 0, 0];
 

  // calculating the length of the midline 
  // nomrlalizing the midline
  let dl = 0;
  for (let i = 0; i < N; i++) {
    let dx = rx[i + 1] - rx[i];
    let dy = ry[i + 1] - ry[i];
    let dz = rz[i + 1] - rz[i];
    dl += Math.sqrt(dx * dx + dy * dy + dz * dz);
  }
  console.log(dl);

 
  // creating tetrahedron vertices
  

  for (var i = 1; i < N +1; i++) {
    

    mid[0] = rx[i];
    mid[1] = ry[i];
    mid[2] = rz[i];

    v[6 * i - 5] = mid[0] - hx[i] * l;
    v[6 * i - 4] = mid[1] - hy[i] * l;
    v[6 * i - 3] = mid[2] - hz[i] * l;
    v[6 * i - 2] = mid[0] + hx[i] * l;
    v[6 * i - 1] = mid[1] + hy[i] * l;
    v[6 * i]     = mid[2] + hz[i] * l;
  }
  for (var i = 0; i < 2*N+1; i++) {
    for (var j = 1; j < 4; j++) {
      v[3*i+j] -= cen[j-1]/n;
      if(j==3){
        v[3*i+j] -=.25;
      }
    };
  };
  // tests
  // console.log(Math.abs(hx[n]*hx[1]+hy[n]*hy[1]+hz[n]*hz[1]+0.925450349781138));
  // tmp = [ hy[n]*hz[1] - hy[1]*hz[n],
  //         hz[n]*hx[1] - hz[1]*hx[n],
  //         hx[n]*hy[1] - hx[1]*hy[n] ];
  // for (var j = 0; j < 3; j++) mid[j] -= tmp[j];
  // console.log([mid[0],mid[1],mid[2]])
  return [v,kappa];
};

var paths = function (n, hl, selector,ind) {
   
  var path = [];
  var v0  = [];
  var v   = [];
  for (var t = 0; t <N+3; t++) {
     v0 = fourierExpansion(n, t, hl)[0];
    v[1] = v0[6 * ind - 5]; 
    v[2] = v0[6 * ind - 4]; 
    v[3] = v0[6 * ind - 3]; 
    v[4] = v0[6 * ind - 2]; 
    v[5] = v0[6 * ind - 1]; 
    v[6] = v0[6 * ind  ];  
    
     if (selector == 1) {
     // path.push(new BABYLON.Vector3(v[1], v[3], -v[2]));
      path.push(
        new BABYLON.Vector3(
          (v[1] + v[4]) / 2,
          (v[3] + v[6]) / 2,
          -(v[2] + v[5]) / 2
        )
      );
    } // corners
    if (selector == 2) {
      path.push(
        new BABYLON.Vector3(
           v[1] ,
           v[3]  ,
          - v[2]  
        )
      );
      // if (i > e / 2 - 1) {
      //   break;
      // }
    } //midpoint
  }
  return path;
};

function poly(N,n, v, i, cm, golfMesh , scheme) {
   var v1 = [v[6 * i - 5], v[6 * i - 3], -v[6 * i - 4]];
  var v2 = [v[6 * i - 2], v[6 * i], -v[6 * i - 1]];
  if (i < N) {
    var v3 = [v[6 * i + 1], v[6 * i + 3], -v[6 * i + 2]];
    var v4 = [v[6 * i + 4], v[6 * i + 6], -v[6 * i + 5]];
  } else {
    var v4 = [v[1], v[3], -v[2]];
    var v3 = [v[4], v[6], -v[5]];
  }

  var pos = v1.concat(v1, v1, v2, v2, v2, v3, v3, v3, v4, v4, v4);

  var indices = [0, 6, 3, 1, 4, 9, 2, 10, 7, 5, 8, 11]; // for mid[j] += tmp[j]/abs;
 
  // var indices = [ 0,3,6, 1,9,4, 2,7,10, 5,11,8 ]; // for mid[j] -= tmp[j]/abs;

  var c1, c2, c3, c4;
   var colors = [];
   if (scheme == 1) {
    // rainbow whole tetrahedron
    c1 =  jet(i); 
     
    
  }   else if (scheme == 2) {
    // Mobius brw 3x
     let n = 3;
      var indx =  Math.floor(i/(Math.round(N/n))) ;
    c1 = jet(Math.round(Math.round(N/n)*indx));
    //  if(indx%2  ==  0){
    //   c1 = [1,0,0,1]
    //   console.log(c1);
    //  }
    //  else if (indx%2  ==  1){
    //    c1 = [0,0,1,1];}

     };
      
      c2  = c1;
     c3 = c1;
     c4 = c1;
     colors = c3.concat(c1, c2, c3, c1, c4, c3, c2, c4, c1, c2, c4);
  
  

  var normals = [];
  BABYLON.VertexData.ComputeNormals(pos, indices, normals);

  var vertexData = new BABYLON.VertexData();

  vertexData.positions = pos;
  vertexData.indices = indices;
  vertexData.colors = colors;
  vertexData.normals = normals;

  vertexData.applyToMesh(cm);
   var ind1 = Math.round(N / (2 * n));

   for (var j=1;j<golfMesh.length+1;j++){
              let l = 2 * j * ind1 - 2 * ind1 + 1;
              var v1 = [v[6 * l - 5], v[6 * l - 3], -v[6 * l - 4]];
              var v2 = [v[6 * l - 2], v[6 * l], -v[6 * l - 1]];
    if(golfMesh[j]){
      golfMesh[j].position  = new BABYLON.Vector3((v1[0] + v2[0]) / 2 ,
      (v1[1] + v2[1]) / 2,(v1[2] + v2[2]) / 2);
    
    };
 
   };

     
}

function anchor(arrow1,tip1, anchorMaterial,v1,v2,hingeLength){ 

 
  // direction for alignment 
  var mod = Math.sqrt((v1[0]-v2[0])**2 + (v1[1]-v2[1])**2 + (v1[2]-v2[2])**2 );
  var direction = new BABYLON.Vector3((v1[0]-v2[0])/mod,(v1[1]-v2[1])/mod,(v1[2]-v2[2])/mod);
  // Compute the rotation to align the body with the given direction
  var axis = BABYLON.Vector3.Cross(BABYLON.Axis.Y, direction);
  var angle = Math.acos(BABYLON.Vector3.Dot(BABYLON.Axis.Y, direction));
  var quater = BABYLON.Quaternion.RotationAxis(axis, angle);
  // position and orientation 
  var mid_position  = new BABYLON.Vector3((v1[0]+v2[0])/2.0,(v1[1]+v2[1])/2.0,(v1[2]+v2[2])/2.0);
  arrow1.position = mid_position;
  arrow1.rotationQuaternion = quater;
  // Assign the height value later
  arrow1.scaling.x =  3.1 * hingeLength;            
  arrow1.scaling.y =  .9 * hingeLength;            
  arrow1.scaling.z =  3.1 * hingeLength;            
  
   
  var h_tip  = .4*hingeLength;
 
  var fac1 = h_tip/2 + hingeLength;
 
 
tip1.rotationQuaternion = quater;
tip1.position  =new BABYLON.Vector3(mid_position.x + fac1*direction.x ,mid_position.y +  fac1*direction.y
                  ,mid_position.z +  fac1* direction.z);
tip1.scaling.x  = h_tip*4.8;
tip1.scaling.y  = h_tip*.8;
tip1.scaling.z  = h_tip*4.8;
  

// tail1.rotationQuaternion = quater;
// tail1.position  =new BABYLON.Vector3(mid_position.x - fac2*direction.x ,mid_position.y -  fac2*direction.y
//                   ,mid_position.z -  fac2*direction.z);
//  tail1.scaling.x  = h_tail*5;
//  tail1.scaling.y  = h_tail;
//  tail1.scaling.z  = h_tail*5;
// Combine the arrow components into a single mesh
//anchor[j] =   BABYLON.Mesh.MergeMeshes([arrow1, tip1, tail1], true, false, null, false, true);

arrow1.material = anchorMaterial;
tip1.material = anchorMaterial;
//tail1.material = anchorMaterial;

};
 


var col = [];
col[1] = [1, 0, 0, 1];
col[2] = [1, 0.25, 0, 1];
col[3] = [1, 0.5, 0, 1];
col[4] = [1, 0.75, 0, 1];
col[5] = [1, 1, 0, 1];
col[6] = [0.75, 1, 0, 1];
col[7] = [0.5, 1, 0, 1];
col[8] = [0.25, 1, 0, 1];
col[9] = [0, 1, 0, 1];
col[10] = [0, 1, 0.25, 1];
col[11] = [0, 1, 0.5, 1];
col[12] = [0, 1, 0.75, 1];
col[13] = [0, 1, 1, 1];
col[14] = [0, 0.75, 1, 1];
col[15] = [0, 0.5, 1, 1];
col[16] = [0, 0.25, 1, 1];
col[17] = [0, 0, 1, 1];
col[18] = [0.25, 0, 1, 1];
col[19] = [0.5, 0, 1, 1];
col[20] = [0.75, 0, 1, 1];
col[21] = [1, 0, 1, 1];
col[22] = [1, 0, 0.75, 1];
col[23] = [1, 0, 0.5, 1];
col[24] = [1, 0, 0.25, 1];

var alpha = 0.1;
var rainbowScheme = function (z) {
  var s = 1 / 6;
  if (z < s) {
    return [1, z / s, 0, alpha];
  } else if (z >= s && z < 2 * s) {
    return [(2 * s - z) / s, 1, 0, alpha];
  } else if (z >= 2 * s && z < 3 * s) {
    return [0, 1, (z - 2 * s) / s, alpha];
  } else if (z >= 3 * s && z < 4 * s) {
    return [0, (4 * s - z) / s, 1, alpha];
  } else if (z >= 4 * s && z < 5 * s) {
    return [(z - 4 * s) / s, 0, 1, alpha];
  } else {
    return [1, 0, (6 * s - z) / s, alpha];
  }
};

var rainbowStepScheme = function (z) {
  var s = 1 / 6;
  if (z < s) {
    return [1, 0, 0, alpha];
  } else if (z >= s && z < 2 * s) {
    return [1, 1, 0, alpha];
  } else if (z >= 2 * s && z < 3 * s) {
    return [0, 1, 0, alpha];
  } else if (z >= 3 * s && z < 4 * s) {
    return [0, 1, 1, alpha];
  } else if (z >= 4 * s && z < 5 * s) {
    return [0, 0, 1, alpha];
  } else {
    return [1, 0, 1, alpha];
  }
};
var brwScheme = function(z) { // black - red - white
  var s = 1/3;
       if (             z <   s ) { return [0,0,0,1]; }
  else if ( z >=   s && z < 2*s ) { return [1,0,0,1]; }
  else                            { return [1,1,1,1]; }
}

var symmetryAxis = function (scene) {
  var d = 0.01;
  var symAxis = BABYLON.Mesh.CreateCylinder(
    "cylinder",
    2,
    d,
    d,
    8,
    1,
    scene,
    false,
    BABYLON.Mesh.DEFAULTSIDE
  );
  var yellow = new BABYLON.StandardMaterial("texture1", scene);
  yellow.diffuseColor = new BABYLON.Color3(1, 1, 0);
  symAxis.material = yellow;
  return symAxis;
};
var meshDispose;
var planeAxis = function (arrowBody ,tip) {

  meshDispose(arrowBody);
   meshDispose(tip);;

  // var l = 1.8;
  // var d = 0.01;
  // var plaAxis = [];
  // plaAxis[1] = BABYLON.Mesh.CreateCylinder(
  //   "cylinder",
  //   l,
  //   d,
  //   d,
  //   8,
  //   1,
  //   scene,
  //   false,
  //   BABYLON.Mesh.DEFAULTSIDE
  // );
  // var yellow = new BABYLON.StandardMaterial("texture1", scene);
  // yellow.diffuseColor = new BABYLON.Color3(1, 1, 0);
  // plaAxis[1].material = yellow;
  // plaAxis[1].rotation.z = Math.PI / 2;
  // for (var i = 2; i < N + 1; ++i) {
  //   plaAxis[i] = plaAxis[1].createInstance("cylinder" + i);
  //   var a = (2 * (i - 1) * Math.PI) / N;
  //   plaAxis[i].rotation.y = a;
  //   plaAxis[i].position.x = (-l / 2) * Math.cos(a);
  //   plaAxis[i].position.z = (l / 2) * Math.sin(a);
  // }
  // plaAxis[1].position.x = -l / 2;
  // return plaAxis;
};

function midlineInit() {
          // Update the path
          var initialPath = [];
          for (var i = 0; i < N + 1; i++) {
            initialPath.push(new BABYLON.Vector3(0, 0, 0));
          }

          optionsMidline = {
            path: initialPath, //vec3 array,
            radius: 0.004, // set the radius of the tube
            updatable: true
          };
          // Update the tube mesh with the new data
          midline = BABYLON.MeshBuilder.CreateTube("midline", optionsMidline, scene);
          midline.material = midlineMaterial;

          // initialize rulings 


          for (var i = 1; i < n + 1; i++) {
            let temp = [];
            temp.push(new BABYLON.Vector3(0, 0, 0));
            temp.push(new BABYLON.Vector3(0, 0, 0));

            optionsRulings = {
              path: temp, //vec3 array,
              radius: 0.01, // set the radius of the tube
              updatable: true
            };

            rulings[i] = BABYLON.MeshBuilder.CreateTube("rulings", optionsRulings, scene);
            rulings[i].material = rulingsMaterial;


          };

        };




// Define a color map based on the index of each point using the jet function
function getJetColor(index, numPoints) {
  const val = index / (numPoints - 1);
  const r = Math.max(
    0,
    Math.min(255, Math.round(255 * (1.5 - 4 * Math.abs(val - 0.5))))
  );
  const g = Math.max(
    0,
    Math.min(
      255,
      Math.round(
        255 * (1.5 - 4 * Math.abs(val - 0.25) - 4 * Math.abs(val - 0.75))
      )
    )
  );
  const b = Math.max(
    0,
    Math.min(255, Math.round(255 * (1.5 - 4 * Math.abs(val - 0.5))))
  );
  return new BABYLON.Color4(r / 255, g / 255, b / 255, 1);
}

function getVibgyorColor(index, numPoints) {
  const step = 1.0 / (numPoints - 1);
  const h = index * step;
  const r = Math.round(
    Math.max(0, Math.min(255, 255 * (1 - Math.abs(h * 6 - 3) - 0.25)))
  );
  const g = Math.round(
    Math.max(0, Math.min(255, 255 * (1 - Math.abs(h * 6 - 2) - 0.25)))
  );
  const b = Math.round(
    Math.max(0, Math.min(255, 255 * (1 - Math.abs(h * 6 - 4) - 0.25)))
  );
  return new BABYLON.Color4(r / 255, g / 255, b / 255, 1);
}


var jet = function(i) {
    
  let r, g, b;

  let N3 = Math.round(N );
  if (i < N3 / 4) {
    r = 0;
    g = Math.floor(4 * i / N3 * 255);
    b = 255;
  } else if (i < N3 / 2) {
    r = 0;
    g = 255;
    b = Math.floor(255 - 4 * (i - N3 / 4) / N3 * 255);
  } else if (i < 3 * N3 / 4) {
    r = Math.floor(4 * (i - N3 / 2) / N3 * 255);
    g = 255;
    b = 0;
  } else {
    r = 255;
    g = Math.floor(255 - 4 * (i - 3 * N3 / 4) / N3 * 255);
    b = 0;
  }


return [  r / 255, g / 255, b / 255, 1]; 
}

var jetK = function(i) {
    
  let r, g, b;

  let N3 = N;
  if (i < N3 / 4) {
    r = 0;
    g = Math.floor(4 * i / N3 * 255);
    b = 255;
  } else if (i < N3 / 2) {
    r = 0;
    g = 255;
    b = Math.floor(255 - 4 * (i - N3 / 4) / N3 * 255);
  } else if (i < 3 * N3 / 4) {
    r = Math.floor(4 * (i - N3 / 2) / N3 * 255);
    g = 255;
    b = 0;
  } else {
    r = 255;
    g = Math.floor(255 - 4 * (i - 3 * N3 / 4) / N3 * 255);
    b = 0;
  }


return `rgba(${r},${g},${b},1)`; 
}



////////////////////////////////
////////// charts /////////////

var renderCanvasK  ;
var contentE ;
var content  ;
var contentK  ;
var contentPath  ;
var contentEWidth;


var x_legend = [];

var legendCount = 13;
var legendSpacing = 40*contentEWidth / (legendCount + 1);

for (var i = 1; i <= legendCount; i++) {
  x_legend.push(i * legendSpacing);
};


var linspace = function (start, end, num) {
  var step = (end - start) / (num - 1);
  var result = [];

  for (var i = 0; i < num; i++) {
    var value = start + step * i;
    result.push(value);
  }

  return result;
}


var getIndexFromLinspace = function (x) {
  let start = xData[0];
  let end = xData[xData.length - 1];
  let length = yData.length;

  if (x < start || x > end || length <= 1) {
    return -1; // x is out of the range or invalid length
  }

  var step = (end - start) / (length - 1);
  var index = Math.round((x - start) / step);

  if (index < 0 || index >= length) {
    return -1; // The calculated index is out of bounds
  }

  return index;
};

var n_to_nui = function (n) {
  let ind = parseInt(Math.round((n - 1) / 2) - 1);
  nu_min = parseInt(parameters[ind][1]);
  nu_opt = parseInt(parameters[ind][2]);
  nu_max = parseInt(parameters[ind][3]);
};

var nui_to_n = function (nu_i) {
  n = null;

  for (let i = 0; i < parameters.length; i++) {
    const row = parameters[i];
    // The second column is at index 1 and the fourth column is at index 3
    if (nu_i >= row[1] && nu_i <= row[3]) {
      n = 2 * i + 3;
      break; // Exit the loop once a matching row is found
    }
  }

  return n; // Will be null if no matching row is found
};



var  deleteDataset = function(chart, label) {
  var datasetIndex = chart.data.datasets.findIndex(dataset => dataset.label === label);

  if (datasetIndex !== -1) {
    chart.data.datasets.splice(datasetIndex, 1);
    chart.update();
  };
};

var  pushDataset = function(chart, ind) {
  var markerCoordinate = { x: xData[ind], y: yData[ind] };


  chart.data.datasets.push({
    label: 'Marker',
    data: [markerCoordinate],
    backgroundColor: 'rgba(0, 0, 0, 0)', // Set the background color to transparent
    borderColor: 'red',
    borderWidth: 3,
    pointRadius: 5,  // Adjust this to change the size of the marker
    pointStyle: 'circle', // optional: to ensure the marker is a circle
    order: 1
  });
  var tempx = [xData[0], xData[ind] + 1.185];
  var tempy = [yData[ind], yData[ind]];
  chart.data.datasets.push({
    label: 'cursor1',
    data: tempx.map((value, index) => ({ x: value, y: tempy[index] })),
    backgroundColor: "rgba(0, 0, 0, 0)",
    borderColor: "rgba(255, 255, 255, 0.5)", // White color with transparency
    borderWidth: 1,
    borderDash: [5, 5], // Dotted line pattern (5px line, 5px space)
    pointRadius: 0,
    order: 2

  });
  var tempx = [xData[ind], xData[ind]];
  var tempy = [35, yData[ind] + 11.5];
  chart.data.datasets.push({
    label: 'cursor2',
    data: tempx.map((value, index) => ({ x: value, y: tempy[index] })),
    backgroundColor: "rgba(0, 0, 0, 0)",
    borderColor: "rgba(255, 255, 255, 0.5)", // White color with transparency
    borderWidth: 1,
    borderDash: [5, 5], // Dotted line pattern (5px line, 5px space)
    pointRadius: 0,
    order: 2

  });

  chart.update();
};


 
var createChart = function(nu_i) {
  if (chart != null) {
    chart.destroy();
  }

  var ctx = document.getElementById('renderCanvasE').getContext('2d', { willReadFrequently: true });

  chart = new Chart(ctx, {
    type: "line",
    data: {
      datasets: [],
    },
    options: {
  
      maintainAspectRatio: false, // Disable aspect ratio to fit the canvas
      scales: {
        x: {
          min: xData[0] - .4,
          max: xData[len_ydata - 1],
          type: "linear",
          position: "bottom",
          ticks: {
            font: {
              size: 16 // Increase the font size for x-axis ticks
            },
            color: 'white',  // this will change the color of the y-axis labels to white
          },
          grid: {
            color: 'rgba(255, 255, 255, 0.1)',  // this will change the color of y-axis grid lines to white (with 10% opacity)
            borderColor: 'white',  // this will change the color of the y-axis itself to white
            //   drawOnChartArea: false,  // Extend gridlines across the entire canvas
            borderDash: [2, 4],  // Optional: Customize the appearance of the gridlines

          },
          title: {
            display: true,
            text: "dimensionless torsion",//'\u03BD', // HTML to display x^2
            color: "white",//
            font: {
              size: 14 // Increase the font size for x-axis ticks
            },
          }
        },
        y: {
          min: 35,
          max: 115,
          type: "linear",
          position: "left",
          ticks: {
            font: {
              size: 16
            },
            color: 'white',
            stepSize: 40,
          },
          grid: {
            color: 'rgba(255, 255, 255, 0.1)',  // this will change the color of y-axis grid lines to white (with 10% opacity)
            borderColor: 'white',  // this will change the color of the y-axis itself to white
            borderDash: [1, 4],  // Optional: Customize the appearance of the gridlines

          },
          title: {
            display: true,
            text: "dimensionless energy", // HTML to display x^2
            color: "white",//
            font: {
              size: 14 // Increase the font size for x-axis ticks
            },
          }
        },

      },
      plugins: {
        legend: {
          display: true,

          labels: {
            filter: function (item, chart) {
              // Return false to hide the legend item
              if (item.text === 'Marker') return false;
              if (item.text === 'cursor1') return false;
              if (item.text === 'cursor2') return false;
              // Return true to show the legend item
              return true;
            },
            color: 'white', // Set the legend labels' color to white
            font: {
              size: 12 // Increase the font size for legend labels
            },
            usePointStyle:true
          },
          //   align: 'start'  // Spread out legend labels to occupy the same space as the horizontal axis

        }
      }
    },
  });
  let ind_marker = nu_i - 1;
  var markerCoordinate = { x: xData[ind_marker], y: yData[ind_marker] };


  chart.data.datasets.push({
    label: 'Marker',
    data: [markerCoordinate],
    backgroundColor: 'rgba(0, 0, 0, 0)', // Set the background color to transparent
    borderColor: 'red',
    borderWidth: 3,
    pointRadius: 5,  // Adjust this to change the size of the marker
    pointStyle: 'circle' // optional: to ensure the marker is a circle
  });

  var tempx = [xData[0], xData[ind_marker] + 1.185];
  var tempy = [yData[ind_marker], yData[ind_marker]];
  chart.data.datasets.push({
    label: 'cursor1',
    data: tempx.map((value, index) => ({ x: value, y: tempy[index] })),
    backgroundColor: "rgba(0, 0, 0, 0)",
    borderColor: "rgba(255, 255, 255, 0.5)", // White color with transparency
    borderWidth: 1,
    borderDash: [5, 5], // Dotted line pattern (5px line, 5px space)
    pointRadius: 0,
    order: 2

  });
  var tempx = [xData[ind_marker], xData[ind_marker]];
  var tempy = [35, yData[ind_marker] + 11.5];
  chart.data.datasets.push({
    label: 'cursor2',
    data: tempx.map((value, index) => ({ x: value, y: tempy[index] })),
    backgroundColor: "rgba(0, 0, 0, 0)",
    borderColor: "rgba(255, 255, 255, 0.5)", // White color with transparency
    borderWidth: 1,
    borderDash: [5, 5], // Dotted line pattern (5px line, 5px space)
    pointRadius: 0,
    order: 2

  });

  var ind = ind_energy.concat();
 // var legendBoxWidths = [20, 20, 20, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300]; // Example: Define desired x-coordinates for each legend label

  for (var i = 0; i < 13; i++) {
    var startIndex = ind[i][1];
    var endIndex = ind[i + 1][1];

    var slicedXData = xData.slice(startIndex, endIndex + 1);
    var slicedYData = yData.slice(startIndex, endIndex + 1);

    var color = `rgb(${colEnergy[i].join(", ")})`;
    chart.data.datasets.push({
      label: "n=" + (2 * i + 3),
      data: slicedXData.map((value, index) => ({ x: value, y: slicedYData[index] })),
      backgroundColor: "rgba(0, 0, 0, 0)",
      borderColor: color,
      borderWidth: 2,
      pointRadius: 0,
      //pointStyle: 'circle' // optional: to ensure the marker is a circle
      order: 2,

      options: {
        // ...
        plugins: {
          title: {
            display: true,
            text: 'curvature of the midline'
          },
          legend: {
            display: true,
            labels: {
              color: 'white',
              font: {
                size: 16
              },
              desiredX: x_legend[i], // Assign the x-location of the legend label from x_legend array
              usePointStyle: true  // Use the same style as points for legend symbols
            },
            align: 'start' // Spread out legend labels to occupy the same space as the horizontal axis
          }
        }
      }
    });
  }
  // Set the boxWidth property for the current legend label
  chart.update(); // Update the chart to display the added datasets
}

function destroyChart(chart) {
  if (chart) {  // Check if chart is defined
    chart.destroy();
    chart = null;
  }
}

// curvature chart

var  createchartK = function(kappa) {
  if (chartK != null) {
    chartK.destroy();
  }
  var kmin = Math.min(...kappa0);
  var kmax = Math.max(...kappa0);

  if (kmax > Math.abs(kmin)) {
    var klim = Math.abs(kmax);
  }
  else {
    var klim = Math.abs(kmin);
  }

  var yticklocation = linspace(-1 * klim, klim, 4);

  var kappa_y = Array.from(kappa).concat();
  //var kappa_x = linspace(1, kappa_y.length, kappa_y.length);
  var kappa_x = linspace(0, 1, kappa_y.length);


  var ctxK = document.getElementById('renderCanvasK').getContext('2d', { willReadFrequently: true });
  chartK = new Chart(ctxK, {
    type: "line",
    data: {
      datasets: [],

    },
    options: {

      responsive: true,
      animation: false, // Disable animation
      maintainAspectRatio: false, // Disable aspect ratio to fit the canvas
      scales: {
        x: {
          min: 0,  
          max: 1,
          type: "linear",
          position: "bottom",
          ticks: {
            color: 'white',
            font: {
              size: 16 // Increase the font size for x-axis ticks
            },  // this will change the color of the y-axis labels to white
          },
          grid: {
            color: 'rgba(255, 255, 255, 0.1)',  // this will change the color of y-axis grid lines to white (with 10% opacity)
            borderColor: 'white',  // this will change the color of the y-axis itself to white
            //   drawOnchartKArea: false,  // Extend gridlines across the entire canvas
            borderDash: [2, 4],  // Optional: Customize the appearance of the gridlines

          },
          title: {
            display: true,
            text: 'arclength', // HTML to display x^2
            color: "white", //
            font: {
              size: 14 // Increase the font size for x-axis ticks
            },
          }
        },
        y: {
          min: -1 * klim - 2,
          max: klim + 2,
          type: "linear",
          position: "left",
          ticks: {
            maxTicksLimit: 4,  // Maximum number of ticks on the y-axis
            color: 'white',  // this will change the color of the y-axis labels to white
            font: {
              size: 16 // Increase the font size for x-axis ticks
            },
            //           callback: function(value, index) {
            //   if (yticklocation.includes(value)) {
            //     return value;
            //   } else {
            //     return '';
            //   }
            // }
          },
          grid: {
            color: 'rgba(255, 255, 255, 0.1)',  // this will change the color of y-axis grid lines to white (with 10% opacity)
            borderColor: 'white',  // this will change the color of the y-axis itself to white
            borderDash: [1, 4],  // Optional: Customize the appearance of the gridlines

          },
          title: {
            display: true,
            text: "curvature of midline", // HTML to display x^2
            color: "white", font: {
              size: 14
            }, //
          }
        },

      },
      plugins: {
        legend: {
          display: false,

          labels: {

            color: 'white', // Set the legend labels' color to white
            font: {
              size: 16 // Increase the font size for legend labels
            },
            usePointStyle:true,
          },

        },

      }
    },
  });
  
  
  // var indices = [];
  // var temp = Math.round(N / (1 * n));
  // for (var i = 1; i < n + 1; i++) {
  //   indices[i - 1] = (i - 1) * temp + Math.round(temp / 2);
  // }
 
 
  // chartK.data.datasets.push({
  //   label: 'markers',
  //   data: indices.map((index) => ({ x: kappa_x[index], y: kappa_y[index] })),
  //   backgroundColor: 'red',
  //   borderColor: 'red',
  //   pointBackgroundColor: 'red',
  //   pointBorderColor: 'red',
  //   pointRadius: 6,
  //   pointHoverRadius: 6,
  //   type: 'scatter',
  //   order: 1,
  // });

  for (var i = 0; i < kappa_y.length - 1; i++) {
    chartK.data.datasets.push({
      label: 'curvature of the midline',
      data: [
        { x: kappa_x[i], y: kappa_y[i] },
        { x: kappa_x[i + 1], y: kappa_y[i + 1] },
      ],
      backgroundColor: 'rgba(0, 0, 0, 0)',
      borderColor: jetK(i),
      borderWidth: 3,
      pointRadius: 0,
      order: 1,
    });
  }

  chartK.update(); // Update the chartK to display the added datasets

}