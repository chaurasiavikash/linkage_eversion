var xmlns='http://www.w3.org/2000/svg'
var s = ' '

// tetrahedron toggle star

function TetTogStar ( id, maxNumOfTets, tetToggled ) {
  
  var star = ebID(id);
  
  var N;
  var nPrevious = 0;
  var nMax = maxNumOfTets;
  
  var state = [];
  for (var i=0; i<nMax+1; ++i) state[i] = 1;
  
  function create( tetv ) {
    var x = [];
    var y = [];
    var angle = 2*Math.PI/N;
    for (var i= 0; i< N+1; ++i) {
      x[i] = Math.sin(i*angle);
      y[i] = Math.cos(i*angle);
    }
    var a = 200;
    var b = 600;
    var p = '';
    for (var i= 0; i< N+1; ++i) {
      var poly = document.createElementNS( xmlns, 'polygon' );
      poly.className.baseVal = 'tetSwitch';
      poly.id = 'tet'+i;
      if ( i == 0 ) { // middle button
        for (var j= 0; j< N+1; ++j) { p = p.concat(a*x[j]+s+a*y[j]+s); }
        poly.setAttribute( 'points', p );
      }
      else {
        poly.setAttribute( 'points', b*x[i-1]+s+b*y[i-1]+s+b*x[i]+s+b*y[i]+s+a*x[i]+s+a*y[i]+s+a*x[i-1]+s+a*y[i-1] );
      }
      star.appendChild(poly);
    }
    for (var i= 0; i< N+1; ++i) {
      var line = document.createElementNS( xmlns, 'polyline' );
      line.className.baseVal = 'separator';
      line.id = 'sep'+i;
      if ( i == 0 ) { // middle button
        line.setAttribute( 'points', p );
        line.style.fill = 'none';
      }
      else {
        line.setAttribute( 'points', a*x[i-1]+s+a*y[i-1]+s+b*x[i-1]+s+b*y[i-1] );
      }
      star.appendChild(line);
    }
    // events /////////////////////////
    for (var i=0; i<N+1; i++) {
      (function(i){
        ebID('tet'+i).addEventListener( 'click', function(){
          if (state[i]==0) {
            if ( i>0 && stateSum()!=N-1 ) {
              ebID('tet'+i).style.fill = '#fff';
              tetv[i].setEnabled(true);
              state[i] = 1;
            }
            else { // all
              for (var j=1; j<N+1; j++) {
                if (state[j]==0) {
                  ebID('tet'+j).style.fill = '#fff';
                  state[j] = 1;
                  tetv[j].setEnabled(true);
                }
              }
              ebID('tet0').style.fill = '#fff';
              state[0] = 1;
            }
          }
          else {
            if (i>0) {
              ebID('tet'+i).style.fill = antra;
              state[i] = 0;
              tetv[i].setEnabled(false);
            }
            else { // all
              for (var j=1; j<N+1; j++) {
                ebID('tet'+j).style.fill = antra;
                state[j] = 0;
                tetv[j].setEnabled(false);
              }
            }
            ebID('tet0').style.fill = antra;
            state[0] = 0;
          }
          document.body.dispatchEvent(tetToggled);
        });
      }(i));
    }
  }
  
  this.changeNumOfTets = function( numOfTets, tetv ) {
    if ( nPrevious != 0 ) {
      for (var i= 0; i<nPrevious+1; ++i) {
        star.removeChild(ebID('tet'+i));
        star.removeChild(ebID('sep'+i));
      }
    }
    N = numOfTets;
    create( tetv );
    if ( stateSum() == 0 ) { // no tetrahedra are shown
      for (var i=0; i<N+1; ++i) {
        ebID('tet'+i).style.fill = antra;
        state[i] = 0;
        if ( i > 0 ) tetv[i].setEnabled(false);
      }
    }
    else for (var i=0; i<N+1; ++i) if (state[i]==0) state[i] = 1; // tetSwitches are white and tets are enabled automatically
    for (var i=N+1; i<nMax+1; ++i) if (state[i]==1) state[i] = 0;
    nPrevious = N;
  }
  
  function stateSum() {
    var tmp = 0;
    for (var i = 0; i<nPrevious+1; ++i) tmp += state[i];
    return tmp;
  }
}
