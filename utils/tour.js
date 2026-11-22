// Automated Shepherd demonstration for the linkage visualizer.
document.addEventListener('DOMContentLoaded', function () {
  const button = document.getElementById('tutorialButton');
  const cursor = document.querySelector('.cursor');
  if (!button || !cursor || !window.Shepherd) return;

  const tour = new Shepherd.Tour({
    tourName: 'linkage-tutorial', exitOnEsc: true, keyboardNavigation: false,
    defaultStepOptions: {
      classes: 'linkage-tutorial', cancelIcon: {enabled: true},
      highlightClass: 'tutorial-target', scrollTo: false, canClickTarget: false,
    },
  });
  window.tour = tour;
  let saved = null, paused = false, closing = false;
  let timer = null, pending = null, generation = 0;
  const api = () => window.linkageTutorial;

  // Draw attention once per tab session, after the visualizer is ready.
  const attentionKey = 'linkage-tutorial-attention-shown';
  let attentionShown = false, attentionTimer = null;
  function stopTutorialAttention() {
    clearTimeout(attentionTimer);
    attentionTimer = null;
    button.classList.remove('tutorial-attention');
  }
  function readyTutorial() {
    button.disabled = false;
    if (attentionShown || saved || closing) return;
    attentionShown = true;
    try {
      if (window.sessionStorage.getItem(attentionKey)) return;
      window.sessionStorage.setItem(attentionKey, '1');
    } catch (_) { /* Storage can be unavailable; keep the first-load highlight. */ }
    button.classList.add('tutorial-attention');
    attentionTimer = setTimeout(stopTutorialAttention, 5000);
  }

  function clearWork() {
    clearTimeout(timer);
    timer = pending = null;
    generation++;
  }

  function later(task, delay) {
    pending = {task, delay};
    if (paused) return;
    timer = setTimeout(() => {
      timer = pending = null;
      task();
    }, delay);
  }

  function buttons(index) {
    return [
      {text: 'Stop', secondary: true, action: () => tour.cancel()},
      {text: paused ? 'Resume' : 'Pause', secondary: true, action: togglePause},
      {text: index === steps.length-1 ? 'Finish' : 'Next', action: () => tour.next()},
    ];
  }

  function togglePause() {
    if (!saved || closing) return;
    paused = !paused;
    clearTimeout(timer);
    timer = null;
    if (!paused && pending) later(pending.task, pending.delay);
    const step = tour.getCurrentStep();
    if (step) step.updateStepOptions({buttons: buttons(tour.steps.indexOf(step))});
  }
  window.toggleTourPause = togglePause;

  const click = id => document.getElementById(id).click();
  const steps = [
    {target: '#buttABCABCABC', title: 'Base sequence 1',
      action: () => click('buttABCABCABC')},
    {target: '#buttABCBCBCAA', title: 'Base sequence 2',
      action: () => click('buttABCBCBCAA')},
    {target: '#buttAAABBBCCC', title: 'Base sequence 3',
      action: () => click('buttAAABBBCCC')},
    {target: '#buttABACBCABC', title: 'Base sequence 4',
      action: () => click('buttABACBCABC')},
    {target: '#sequenceInput', title: 'Custom sequence',
      text: 'Enter any order containing exactly three A, three B, and three C. Here is an example.',
      action: () => {
        const input = document.getElementById('sequenceInput');
        input.value = 'CBACBACBA';
        input.dispatchEvent(new Event('input', {bubbles: true}));
      }},
    {target: '#generateMotion', title: 'Generate motion',
      text: 'Generate the motion for your sequence directly in the browser.',
      action: async () => {
        document.getElementById('sequenceForm').requestSubmit();
        await api().whenIdle();
      }},
    {target: '#buttSwitch', title: 'Swap adjacent links',
      text: 'Each click exchanges the next neighboring pair and updates the motion.',
      action: () => click('buttSwitch')},
    {target: '#buttPlayPause', title: 'Play and pause',
      text: 'Play the motion or pause it. You can also drag the surrounding clock ring to choose a phase.',
      action: () => api().setPlaying(true)},
    {target: '#speedSlider', title: 'Playback speed',
      text: 'Adjust playback from the default speed to twice that speed.',
      action: () => api().setSpeed(200)},
    {target: '#buttFullScreen', title: 'Full screen',
      text: 'Use Full screen when you want more space for the visualizer.'},
    {target: '#buttShowMidVol', title: 'Swept volume',
      text: 'Show the surface swept out by the linkage as it moves.',
      action: () => api().setSweptVolume(true)},
    {target: '#buttShowMid', title: 'Hinge trajectories',
      text: 'Yellow curves show the paths traced by hinge centers during eversion.',
      action: () => { api().setSweptVolume(false); api().setTrajectories(true); }},
    {target: '#exportData', title: 'Export data',
      text: 'After the tutorial, use Export data to pause and save the configuration, full motion, trajectories, and Python plotting scripts.',
      action: () => api().setPlaying(false)},
    {target: '#leftButton', title: 'Left-handed linkage',
      text: 'Switch chirality to inspect the left-handed linkage.',
      action: () => click('leftButton')},
    {target: '#rightButton', title: 'Right-handed linkage',
      text: 'Switch back to the right-handed linkage.',
      action: () => click('rightButton')},
    {target: '#buttToggleView', title: 'View along the axis',
      text: 'Use this button to alternate between the two aligned views.',
      action: () => document.getElementById('buttToggleView').dispatchEvent(new CustomEvent('alignView', {detail: 'axis'}))},
    {target: '#buttToggleView', title: 'Second aligned view',
      text: 'The next click shows the other view of the linkage.',
      action: () => document.getElementById('buttToggleView').dispatchEvent(new CustomEvent('alignView', {detail: 'plane'}))},
    {target: '#renderCanvas', attach: '#buttToggleView', title: 'Explore in 3D',
      text: 'Drag the canvas to orbit the linkage. Scroll or pinch to zoom.',
      action: () => api().orbit()},
    {target: '#backgroundButton', title: 'Background color',
      text: 'Choose a background color that suits your display or presentation. Finish restores your starting motion and view.'},
  ];

  steps.forEach((demo, index) => {
    const target = document.querySelector(demo.target);
    const side = target.closest('.leftbar') ? 'right' : target.closest('.rightbar') ? 'left' : 'top';
    tour.addStep({
      id: 'linkage-step-'+index,
      title: demo.title,
      text: demo.text,
      attachTo: {element: demo.attach || demo.target, on: demo.attach ? 'left' : side},
      buttons: buttons(index),
      when: {
        show() {
          clearWork();
          const token = generation;
          this.updateStepOptions({buttons: buttons(index)});
          target.scrollIntoView({block: 'nearest', inline: 'nearest'});
          const rect = target.getBoundingClientRect();
          cursor.classList.remove('expand');
          cursor.style.left = (rect.left+rect.width/2)+'px';
          cursor.style.top = (rect.top+rect.height/2)+'px';
          later(async () => {
            try {
              if (demo.action) {
                cursor.classList.add('expand');
                await demo.action();
              }
              if (saved && token === generation) later(() => tour.next(), 3600);
            } catch (error) {
              console.error('Tutorial step failed:', error);
              tour.cancel();
            }
          }, 850);
        },
        hide: clearWork,
      },
    });
  });

  async function finish() {
    clearWork();
    if (!saved || closing) return;
    closing = true;
    const initial = saved;
    saved = null;
    button.disabled = true;
    document.body.classList.remove('tutorial-running');
    cursor.classList.remove('expand');
    cursor.style.cssText = initial.cursorStyle;
    try {
      await api().restore(initial.motion);
      initial.scroll.forEach(([element, top]) => { element.scrollTop = top; });
    } finally {
      paused = closing = false;
      button.disabled = false;
      button.textContent = 'Tutorial';
      button.setAttribute('aria-pressed', 'false');
    }
  }
  tour.on('cancel', finish);
  tour.on('complete', finish);

  window.startTour = function () {
    if (!api() || saved || closing) return;
    stopTutorialAttention();
    saved = {
      motion: api().capture(), cursorStyle: cursor.style.cssText,
      scroll: Array.from(document.querySelectorAll('.leftbar, .rightbar'), element => [element, element.scrollTop]),
    };
    paused = false;
    api().setPlaying(false);
    api().setSweptVolume(false);
    api().setTrajectories(false);
    document.body.classList.add('tutorial-running');
    button.textContent = 'Stop tutorial';
    button.setAttribute('aria-pressed', 'true');
    tour.start();
  };
  button.addEventListener('click', () => saved ? tour.cancel() : window.startTour());
  document.addEventListener('linkageTutorialReady', readyTutorial);
  if (api()) readyTutorial();
});
