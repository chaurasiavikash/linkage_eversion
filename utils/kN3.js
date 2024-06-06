var nu_min  ;
var nu_max  ;
var nu_opt  ;

var numData;
var parameters;
var nu_i;
var N;
var N1;

var Pi = Math.PI;
 var n;
var h;
var kappa0 = [];  // curvature of the midline at time t = 0
var kappa = [];  // curvature of the midline at time t = 0

var plotClicked;
var bx_data;
var by_data;
var bz_data;
var tmax;



var init_parameters = function (sequence) {
  var strData;
  if(sequence==1){
    strData = 'ACBACBACB';
  //  tmax = 271;
  } 
  else if(sequence==2){
    strData = 'ABCABCABC';
    //tmax = 270;
  }
  else if(sequence==3){
    strData = 'AAABCBCBC';
    //tmax =   255;
  }
  else{
    strData = 'AAACCCBBB';
    //tmax =   254;
  }
   
  bx_data = numData[strData].bx_data  ;
  by_data = numData[strData].by_data  ;
  bz_data = numData[strData].bz_data  ;
   tmax = bx_data[1].length;
  // creating the init_parameters
};
var avg = function (arr) {
   return arr.reduce((a, b) => a + b, 0) / arr.length;
};
  

var cen = [0, 0, 0];
 

//  midline and binormal loaded from the saved file 
const arrayColumn = (arr, n) => arr.map(x => x[n]);

//the function below generates animation data 
var animationData = function(j) {
  let ind = j % (tmax); // ensuring that index is not out of bound

  hx = arrayColumn(bx_data,ind); //
  hy = arrayColumn(by_data,ind); //
  hz = arrayColumn(bz_data,ind); //
  
  let tx = [];
  let ty = [];
  let tz = []; 
 
  let h = 1;
     for (var i = 1; i < N + 1; i++) {
    

    tx[i] = (hy[i - 1] * hz[i] - hy[i] * hz[i - 1])/tau/h;
    ty[i] = (hz[i - 1] * hx[i] - hz[i] * hx[i - 1])/tau/h;
    tz[i] = (hx[i - 1] * hy[i] - hx[i] * hy[i - 1])/tau/h;
 

  }
  // if (isNaN(hx  ) || isNaN(hy  ) || isNaN(hz   )) {
  //   console.log(`Found NaN at index ${i}`);
  // }

  tx[0] = tx[N];
  ty[0] = ty[N];
  tz[0] = tz[N];
 
///////////////////////


 
  // Midline
  rx[0] = 0.0;
  ry[0] = 0.0;
  rz[0] = 0.0;
  h = 1;
  for (var i = 0; i < N; i++) {
    rx[i + 1] = rx[i] + h * tx[i + 1];
    ry[i + 1] = ry[i] + h * ty[i + 1];
    rz[i + 1] = rz[i] + h * tz[i + 1];
  }

  cen[0] = avg(rx);
  cen[1] = avg(ry);
  cen[2] = avg(rz);

  for (var i = 0; i < N + 1; i++) {
    rx[i] = rx[i] - cen[0];
    ry[i] = ry[i] - cen[1];
    rz[i] = rz[i] - cen[2];
  }
  // nomrlalizing the midline
  let dl = 0;
  for (let i = 0; i < N; i++) {
    let dx = rx[i + 1] - rx[i];
    let dy = ry[i + 1] - ry[i];
    let dz = rz[i + 1] - rz[i];
    dl += Math.sqrt(dx * dx + dy * dy + dz * dz);
  }
  dl = dl / 6;
  for (let i = 0; i <= N; i++) {
    rx[i] /= dl;
    ry[i] /= dl;
    rz[i] /= dl;
  } 

  return [hx,hy,hz,rx,ry,rz];
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

   kappa = rx;
  
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

  let len_path = Math.round(tmax/2.0)+2;
   for (var t = 1; t < len_path; t++) {
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
 
// assigning position and orientation to arrows, connectors etc. 

function locationAndDirection(v1, v2) {
  var mod = Math.sqrt((v1[0] - v2[0])**2 + (v1[1] - v2[1])**2 + (v1[2] - v2[2])**2);
  let direction = new BABYLON.Vector3((v1[0] - v2[0]) / mod, (v1[1] - v2[1]) / mod, (v1[2] - v2[2]) / mod);
  let mid_position = new BABYLON.Vector3((v1[0] + v2[0]) / 2.0, (v1[1] + v2[1]) / 2.0, (v1[2] + v2[2]) / 2.0);
  let lenV1V2 = mod;
 
  return [  mid_position, direction];
}

function assignStructure(v,hingeLength,arrowBody,tip,topDisc,bottomDisc,hingeBody,connector){
  
  
  for (var i = 1;i<N+1;++i){

    var v1 = [v[6 * i - 5], v[6 * i - 3], -v[6 * i - 4]];
    var v2 = [v[6 * i - 2], v[6 * i], -v[6 * i - 1]];
 
   // assigning coordinate and direction to arrows 
   // direction for alignment  
  // Compute the rotation to align the body with the given direction

  let lenV1V2 = Math.sqrt((v1[0] - v2[0])**2 + (v1[1] - v2[1])**2 + (v1[2] - v2[2])**2);

  let [  mid_position, direction] = locationAndDirection(v1, v2);

 
  var axis = BABYLON.Vector3.Cross(BABYLON.Axis.Y, direction);
  var angle = Math.acos(BABYLON.Vector3.Dot(BABYLON.Axis.Y, direction));
  var quater = BABYLON.Quaternion.RotationAxis(axis, angle);
  // position and orientation 
 
  arrowBody[i].position = mid_position;
  arrowBody[i].rotationQuaternion = quater;
  // Assign the height value later
  arrowBody[i].scaling.x =  3.1 * hingeLength;            
  arrowBody[i].scaling.y =  .8 * hingeLength;            
  arrowBody[i].scaling.z =  3.1 * hingeLength;            
  
   
  var h_tip  = .5*hingeLength;
 
  var fac1 = h_tip/2 + .9*hingeLength;
 
 
tip[i].rotationQuaternion = quater;
tip[i].position  =new BABYLON.Vector3(mid_position.x + fac1*direction.x ,mid_position.y +  fac1*direction.y
                  ,mid_position.z +  fac1* direction.z);
tip[i].scaling.x  = h_tip*4.8;
tip[i].scaling.y  = h_tip*.8;
tip[i].scaling.z  = h_tip*4.8;
  
 // assigning hinges 

 
 hingeBody[i].position = mid_position;
 hingeBody[i].rotationQuaternion = quater;
 hingeBody[i].scaling.x =  4.1 * hingeLength;
 hingeBody[i].scaling.y =  .4 * hingeLength;
 hingeBody[i].scaling.z =  4.1 * hingeLength;
 
 fac1 = hingeLength;

 topDisc[i].rotationQuaternion = quater;
 topDisc[i].position  =new BABYLON.Vector3(mid_position.x + fac1*direction.x ,mid_position.y +  fac1*direction.y
                   ,mid_position.z +  fac1* direction.z);
 topDisc[i].scaling.x  = 5.1 * hingeLength;
 topDisc[i].scaling.y  = 0;
 topDisc[i].scaling.z  = 5.1 * hingeLength;

 fac1 = -hingeLength;

 bottomDisc[i].rotationQuaternion = quater;
 bottomDisc[i].position  =new BABYLON.Vector3(mid_position.x + fac1*direction.x ,mid_position.y +  fac1*direction.y
                   ,mid_position.z +  fac1* direction.z);
 bottomDisc[i].scaling.x  = 5.1 * hingeLength;
 bottomDisc[i].scaling.y  = 0;
 bottomDisc[i].scaling.z  = 5.1 * hingeLength;


// assigning connectors 
var iPlus1;
if(i<N){
  iPlus1 = i+1;
   var v3 = [v[6 * iPlus1 - 5], v[6 * iPlus1 - 3], -v[6 * iPlus1 - 4]];
  var v4 = [v[6 * iPlus1 - 2], v[6 * iPlus1], -v[6 * iPlus1 - 1]];
}
else{
  iPlus1 = 1;

  var v4 = [v[6 * iPlus1 - 5], v[6 * iPlus1 - 3], -v[6 * iPlus1 - 4]];
  var v3 = [v[6 * iPlus1 - 2], v[6 * iPlus1], -v[6 * iPlus1 - 1]];
};
 
   // assigning coordinate and direction to arrows 
   // direction for alignment  
  // Compute the rotation to align the body with the given direction

  let temp1 = [(v1[0] + v2[0]) / 2, (v1[1] + v2[1]) / 2, (v1[2] + v2[2]) / 2];
  let temp2 = [(v3[0] + v4[0]) / 2, (v3[1] + v4[1]) / 2, (v3[2] + v4[2]) / 2];
    
   [mid_position, direction] = locationAndDirection(temp1,temp2);

  axis = BABYLON.Vector3.Cross(BABYLON.Axis.Y, direction);
  angle = Math.acos(BABYLON.Vector3.Dot(BABYLON.Axis.Y, direction));
  quater = BABYLON.Quaternion.RotationAxis(axis, angle);


  connector[i].position = mid_position;
  connector[i].rotationQuaternion = quater;
   
  x1 = (v1[0]+v2[0])/2;
  y1 = (v1[1]+v2[1])/2;
  z1 = (v1[2]+v2[2])/2;

  x2 = (v3[0]+v4[0])/2;
  y2 = (v3[1]+v4[1])/2;
  z2 = (v3[2]+v4[2])/2;

   lenV1V2 = 0.5*Math.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2);

   connector[i].scaling.x =  .2;            
   connector[i].scaling.y =  1*lenV1V2;;            
   connector[i].scaling.z =  .2;            
   

  };

  /// assigning coordinates to the N+1th arrow 
  i = 1;
  var v2 = [v[6 * i - 5], v[6 * i - 3], -v[6 * i - 4]];
    var v1 = [v[6 * i - 2], v[6 * i], -v[6 * i - 1]];
 
    i = N+1;
   // assigning coordinate and direction to arrows 
   // direction for alignment  
  // Compute the rotation to align the body with the given direction

  let [mid_position,direction] = locationAndDirection(v1,v2) ;
 
  var axis = BABYLON.Vector3.Cross(BABYLON.Axis.Y, direction);
  var angle = Math.acos(BABYLON.Vector3.Dot(BABYLON.Axis.Y, direction));
  var quater = BABYLON.Quaternion.RotationAxis(axis, angle);
  // position and orientation 
 
  arrowBody[i].position = mid_position;
  arrowBody[i].rotationQuaternion = quater;
  // Assign the height value later
  arrowBody[i].scaling.x =  3.1 * hingeLength;            
  arrowBody[i].scaling.y =  .8 * hingeLength;            
  arrowBody[i].scaling.z =  3.1 * hingeLength;            
  
   
  var h_tip  = .5*hingeLength;
 
  var fac1 = h_tip/2 + .9*hingeLength;
 
 
tip[i].rotationQuaternion = quater;
tip[i].position  =new BABYLON.Vector3(mid_position.x + fac1*direction.x ,mid_position.y +  fac1*direction.y
                  ,mid_position.z +  fac1* direction.z);
tip[i].scaling.x  = h_tip*4.8;
tip[i].scaling.y  = h_tip*.8;
tip[i].scaling.z  = h_tip*4.8;

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
 