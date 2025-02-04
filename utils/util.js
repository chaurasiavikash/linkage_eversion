var antra = '#555';


var hx ;
var hy;
var hz;

var rx;
var ry;
var rz;


function FormatNumberLength(num, length) {
    var r = "" + num;
    while (r.length < length) {
        r = "\xa0" + r;
    }
    return r;
}

function buttonSwitchColor(butt,on) {
  if ( on ) {
    butt.style.color = antra;
    butt.style.backgroundColor = '#fff';
  }
  else {
    butt.style.color = '#fff';
    butt.style.backgroundColor = antra;
  }
}

function sum(v) {
  var tmp = 0;
  for (var i = 0; i < v.length; i++) {
    tmp += v[i];
  }
  return tmp;
}

var pi = Math.PI;
var sqr2 = Math.sqrt(2);
var sqr3 = Math.sqrt(3);
var sqr6 = Math.sqrt(6);
function sin(x) { return Math.sin(x) }
function cos(x) { return Math.cos(x) }
function tan(x) { return Math.tan(x) }

function tripod(scene,scale,translate,alpha)
{
  var l = 0.15*scale;
  var d = 0.015*scale;
  var cyl1 = BABYLON.Mesh.CreateCylinder('cylinder',l,d,d,4,1, scene, false, BABYLON.Mesh.DEFAULTSIDE);
  var cyl2 = BABYLON.Mesh.CreateCylinder('cylinder',l,d,d,4,1, scene, false, BABYLON.Mesh.DEFAULTSIDE);
  var cyl3 = BABYLON.Mesh.CreateCylinder('cylinder',l,d,d,4,1, scene, false, BABYLON.Mesh.DEFAULTSIDE);
  var red   = new BABYLON.StandardMaterial('texture1', scene);
  var green = new BABYLON.StandardMaterial('texture1', scene);
  var blue  = new BABYLON.StandardMaterial('texture1', scene);
  red.diffuseColor    = new BABYLON.Color3(1,0,0);
  green.diffuseColor  = new BABYLON.Color3(0,1,0);
  blue.diffuseColor   = new BABYLON.Color3(0,0.4,1);
  if ( alpha ) {
    red.alpha = alpha;
    green.alpha = alpha;
    blue.alpha = alpha;
  }
  cyl1.material = red;
  cyl1.rotation.z = - Math.PI/2;
  cyl1.position.x = l/2;
  cyl2.material = green;
  cyl2.rotation.x = - Math.PI/2;
  cyl2.position.z = - l/2;
  cyl3.material = blue;
  cyl3.position.y = l/2;
  if ( translate ) {
    cyl1.position.x += translate[0]; cyl2.position.x += translate[0]; cyl3.position.x += translate[0];
    cyl1.position.y += translate[1]; cyl2.position.y += translate[1]; cyl3.position.y += translate[1];
    cyl1.position.z += translate[2]; cyl2.position.z += translate[2]; cyl3.position.z += translate[2];
  }
};

var ebID = function(id) {
  return document.getElementById(id);
}

// custom slider //////////////////////////////////////////////////////////////
var sliderEvent;
sliderEvent = document.createEvent('HTMLEvents');
sliderEvent.initEvent('sliderChange', true, true);
sliderEvent.eventName = 'sliderChange';

function Slider ( id, min, max, init, prefix, digits, suffix ) {
  var slider = ebID(id);
  var val;
  var maxi = max;
  var pt = slider.createSVGPoint();
  function cursorPoint(ev){
    
    return 
  }
  var dragging = false;
  slider.addEventListener('mousedown',startDragMouse,false);
  slider.addEventListener('touchstart',startDragTouch,false);
  document.body.addEventListener('mousemove',dragMouse,false);
  document.body.addEventListener('touchmove',dragTouch,false);
  window.addEventListener('mouseup',stopDrag,false);
  window.addEventListener('touchend',stopDrag,false);
  function startDragMouse(ev) {
    dragging = true;
    dragMouse(ev);
    slider.style.opacity = 1;
  }
  function startDragTouch(ev) {
    dragging = true;
    dragTouch(ev);
    slider.style.opacity = 1;
  }
  function dragMouse(ev) {
    if(dragging) {
      pt.x = ev.clientX;
      pt.y = ev.clientY;
      var loc = pt.matrixTransform(slider.getScreenCTM().inverse());
      var x = Math.min(Math.max(loc.x,-80),80);
      update(x);
      slider.dispatchEvent(sliderEvent);
    }
  }
  function dragTouch(ev) {
    if(dragging) {
      pt.x = ev.touches[0].clientX;
      pt.y = ev.touches[0].clientY;
      var loc = pt.matrixTransform(slider.getScreenCTM().inverse());
      var x = Math.min(Math.max(loc.x,-80),80);
      update(x);
      slider.dispatchEvent(sliderEvent);
    }
  }
  function stopDrag() {
    dragging = false;
    slider.style.opacity = '';
  }
  function update(x) {
    val = min + (x+80)/160*(maxi-min);
    
    if(id=='nuSlider'){
    val = Math.round(val); // Round val to the nearest integer
    };
    slider.getElementsByTagName('rect')[0].setAttribute('x',x-19);
    slider.getElementsByTagName('text')[0].innerHTML = prefix+val.toFixed(digits)+suffix;

  }
  update((init-min)/(maxi-min)*160-80); // initial
  this.getValue = function() { return val; }
  this.setValue = function(va) { val = va; }
  this.setMax = function(max) { maxi = max; }
}
///////////////////////////////////////////////////////////////////////////////
var t;
function test(scene)
{
  var s1 = BABYLON.Mesh.CreateSphere('sphere1', 32, 0.1, scene);
 // var s2 = BABYLON.Mesh.CreateSphere('sphere1', 32, 0.1, scene);
  //var s3 = BABYLON.Mesh.CreateSphere('sphere1', 32, 0.1, scene);
  var red   = new BABYLON.StandardMaterial('texture1', scene);
  var green = new BABYLON.StandardMaterial('texture1', scene);
  var blue  = new BABYLON.StandardMaterial('texture1', scene);
  red.diffuseColor    = new BABYLON.Color3(1,0,0);
  green.diffuseColor  = new BABYLON.Color3(0,1,0);
  blue.diffuseColor   = new BABYLON.Color3(0,0.3,1);
  s1.material = green;

  
   
 

  s1.position.x =rx[0];
  s1.position.y =ry[0];
  s1.position.z =rz[0];
  // s2.material = green;
  // s2.position.y = 0;
  // s3.material = blue;
  // s3.position.z = 0;
};

///////////////////////////////////////////////////////////////////////////////

function toggleFullScreen() {
  var element = document.body;
  if (!document.fullscreenElement && // if not full screen
      !document.mozFullScreenElement &&
      !document.webkitFullscreenElement &&
      !document.msFullscreenElement ) {
    var requestMethod = element.requestFullScreen || element.webkitRequestFullScreen
                  || element.mozRequestFullScreen || element.msRequestFullScreen;
    if (requestMethod) requestMethod.call(element);
  }
  else {
         if (document.exitFullscreen      ) document.exitFullscreen();
    else if (document.msExitFullscreen    ) document.msExitFullscreen();
    else if (document.mozCancelFullScreen ) document.mozCancelFullScreen();
    else if (document.webkitExitFullscreen) document.webkitExitFullscreen();
  }
}

// remove hover functionality if a touch event is registered
var css = '.button:hover,.touchElement:hover,.tetSwitch:hover,.slider:hover {opacity:0.7;} #buttNumOfTetMinus:hover ~ #infoTet {opacity:0.7;} #buttNumOfTetPlus:hover ~ #infoTet {opacity:0.7;}';
var style = document.createElement('style');
document.addEventListener('touchstart', function(event) {
  if (style.styleSheet) {
    style.styleSheet.cssText = css;
  } else {
    style.appendChild(document.createTextNode(css));
  }
  document.getElementsByTagName('head')[0].appendChild(style);
} );

// use esc-key to toggle full screen mode
document.addEventListener('keydown', function(event) { if ( event.keyCode == '27' ) toggleFullScreen(); } );
