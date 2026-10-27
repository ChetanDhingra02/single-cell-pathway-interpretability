(function () {
  'use strict';
  const D = window.__DATA;
  const $ = (s, r) => (r || document).querySelector(s);
  const $$ = (s, r) => Array.from((r || document).querySelectorAll(s));
  const body = document.body;
  const reduce = window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  const comma = n => Number(n).toLocaleString('en-US');
  const f3 = v => (v === null || v === undefined) ? '–' : Number(v).toFixed(3);
  const pct = v => Math.round(v * 100) + '%';
  const sgn = v => (v >= 0 ? '+' : '−') + Math.abs(v).toFixed(3);
  const clamp = (v, a, b) => Math.max(a, Math.min(b, v));
  const ease = t => 1 - Math.pow(1 - t, 3);
  const easeIO = t => t < .5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2;
  const isTech = () => body.classList.contains('tech');

  const UNITS = D.units.slice();
  const byId = Object.fromEntries(UNITS.map(u => [u.id, u]));
  const IFNA = 'HALLMARK_INTERFERON_ALPHA_RESPONSE';
  const least = UNITS.slice().sort((a, b) => a.p - b.p)[0];
  const purest = UNITS.slice().sort((a, b) => b.p - a.p)[0];

  const CH = [
    { id: 'question', n: '01', title: 'The question', teaser: 'What does a pathway label really guarantee?', fill: '#dde8f4', line: '#3d6fb6' },
    { id: 'method', n: '02', title: 'How we tested', teaser: 'Three faithfulness checks, then external validation', fill: '#e8e3f5', line: '#7e6bc4' },
    { id: 'faith', n: '03', title: 'Partly faithful', teaser: 'Median attribution purity: ' + pct(D.faith.median_purity), fill: '#f8e5d6', line: '#c8693a' },
    { id: 'overlap', n: '04', title: 'Pathway overlap', teaser: 'More overlap is linked to more correlated units', fill: '#fbf0c7', line: '#a88418' },
    { id: 'fix', n: '05', title: 'The split', teaser: 'Separate shared genes from pathway-specific genes', fill: '#e3ecdf', line: '#3f7d58' },
    { id: 'external', n: '06', title: 'External validation', teaser: 'Independent validation on five external donors', fill: '#dff0e5', line: '#3f8f7a' },
    { id: 'takeaway', n: '07', title: 'Takeaway', teaser: 'What the results actually support', fill: '#f6dfdf', line: '#b5525b' },
  ];

  /* ---------------- tooltip ---------------- */
  const tip = $('#tip');
  function showTip(html, ev) {
    tip.innerHTML = html; tip.classList.add('on');
    const app = $('#app').getBoundingClientRect();
    let x = ev.clientX - app.left + 14, y = ev.clientY - app.top - tip.offsetHeight - 12;
    if (x + tip.offsetWidth > app.width - 8) x = ev.clientX - app.left - tip.offsetWidth - 14;
    if (y < 8) y = ev.clientY - app.top + 18;
    tip.style.left = x + 'px'; tip.style.top = y + 'px';
  }
  const hideTip = () => tip.classList.remove('on');

  /* ---------------- counters ---------------- */
  function countUp(el, to, dec, dur) {
    dec = dec || 0; dur = dur || 1200;
    if (reduce) { el.textContent = to.toFixed(dec); return; }
    const from = parseFloat(el.dataset.cur || '0'); const t0 = performance.now();
    (function step(t) {
      const k = Math.min(1, (t - t0) / dur), v = from + (to - from) * ease(k);
      el.textContent = v.toFixed(dec);
      if (k < 1) requestAnimationFrame(step); else el.dataset.cur = String(to);
    })(t0);
  }

  /* ---------------- static facts ---------------- */
  $$('[data-fact]').forEach(e => { e.textContent = comma(D.facts[e.dataset.fact]); });
  $$('[data-fact-r]').forEach(e => { e.textContent = D.recon[e.dataset.factR].toFixed(4); });
  $('#faith-pur').textContent = D.faith.median_purity.toFixed(3);
  $('#faith-below').textContent = D.faith.n_below_half;
  $('#faith-sp').textContent = D.faith.median_spearman_own.toFixed(3);
  $('#faith-below-pct').textContent = pct(D.faith.share_below_half);
  $('#faith-pct-txt').textContent = pct(D.faith.median_purity);
  $$('[data-ms]').forEach(e => { e.textContent = D.modules[e.dataset.ms]; });
  $$('[data-ov]').forEach(e => { const [m, c, k] = e.dataset.ov.split('|'); e.textContent = f3(D.overlap.summary[m + '|' + c][k]); });
  $$('[data-vd]').forEach(e => { const v = D.verdicts[e.dataset.vd]; e.textContent = sgn(v.med) + ' (' + v.nb + '/' + v.n + ')'; });
  $('#faith-abl').textContent = Math.round(D.faith.median_ablation_ratio);
  $('#faith-abl-t').textContent = D.faith.median_ablation_ratio.toFixed(1);
  $('#ctl-d').textContent = f3(D.control.d_own.med); $('#ctl-e').textContent = f3(D.control.e_own.med); $('#ctl-o').textContent = f3(D.control.e_orig.med);

  /* ================================================================
     HOME: the 47 pathway dials as living cells (cursor pushes them)
     ================================================================ */
  const home = (function () {
    const cv = $('#home-cells'), dish = $('#home-dish'), ctx = cv.getContext('2d');
    let W = 0, R = 0, cells = [], mouse = null, running = false, hover = -1, t0 = performance.now();
    function size() {
      const dpr = window.devicePixelRatio || 1; W = dish.clientWidth; R = W / 2;
      cv.width = W * dpr; cv.height = W * dpr; ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      const r = Math.max(7, W * 0.026);
      if (!cells.length) {
        cells = UNITS.map((u, i) => {
          const a = Math.random() * 6.283, d = Math.sqrt(Math.random()) * (R * 0.72);
          return { u, x: R + Math.cos(a) * d, y: R + Math.sin(a) * d, vx: (Math.random() - .5) * .6, vy: (Math.random() - .5) * .6, ph: Math.random() * 6.28, r };
        });
      } else cells.forEach(c => { c.r = r; c.x = clamp(c.x, r, W - r); c.y = clamp(c.y, r, W - r); });
    }
    function step() {
      const now = performance.now(), t = (now - t0) / 1000, inner = R - 16;
      for (const c of cells) {
        c.vx += (Math.random() - .5) * .05; c.vy += (Math.random() - .5) * .05;
        if (mouse) {
          const dx = c.x - mouse.x, dy = c.y - mouse.y, d = Math.hypot(dx, dy), reach = W * 0.16;
          if (d < reach && d > 0.1) { const f = (1 - d / reach) * 1.1; c.vx += dx / d * f; c.vy += dy / d * f; }
        }
        c.vx *= .94; c.vy *= .94;
        c.x += c.vx; c.y += c.vy;
        const dx = c.x - R, dy = c.y - R, d = Math.hypot(dx, dy), lim = inner - c.r;
        if (d > lim) { c.x = R + dx / d * lim; c.y = R + dy / d * lim; c.vx -= dx / d * .6; c.vy -= dy / d * .6; }
      }
      for (let i = 0; i < cells.length; i++) for (let j = i + 1; j < cells.length; j++) {
        const a = cells[i], b = cells[j], dx = b.x - a.x, dy = b.y - a.y, d = Math.hypot(dx, dy), m = a.r + b.r + 2;
        if (d < m && d > 0.01) { const p = (m - d) / 2 / d; a.x -= dx * p; a.y -= dy * p; b.x += dx * p; b.y += dy * p; }
      }
      ctx.clearRect(0, 0, W, W);
      cells.forEach((c, i) => {
        const low = c.u.p < 0.5, w = Math.sin(t * 1.6 + c.ph) * 0.07, rr = i === hover ? c.r * 1.35 : c.r;
        ctx.beginPath(); ctx.ellipse(c.x, c.y, rr * (1 + w), rr * (1 - w), c.ph + t * .2, 0, 6.283);
        ctx.fillStyle = low ? 'rgba(243,205,179,.92)' : 'rgba(196,215,239,.95)'; ctx.fill();
        ctx.lineWidth = i === hover ? 2.2 : 1.4; ctx.strokeStyle = low ? '#c8693a' : '#3d6fb6'; ctx.stroke();
        ctx.beginPath(); ctx.arc(c.x + rr * .2, c.y - rr * .15, rr * .36, 0, 6.283);
        ctx.fillStyle = low ? 'rgba(200,105,58,.6)' : 'rgba(61,111,182,.55)'; ctx.fill();
      });
      if (running) requestAnimationFrame(step);
    }
    dish.addEventListener('mousemove', e => {
      const r = dish.getBoundingClientRect(); mouse = { x: e.clientX - r.left, y: e.clientY - r.top };
      let best = -1, bd = 1e9;
      cells.forEach((c, i) => { const d = Math.hypot(c.x - mouse.x, c.y - mouse.y); if (d < c.r + 6 && d < bd) { bd = d; best = i; } });
      hover = best;
      if (best >= 0) showTip('<b>' + cells[best].u.n + '</b>' + pct(cells[best].u.p) + ' of its attribution on its own pathway genes', e); else hideTip();
    });
    dish.addEventListener('mouseleave', () => { mouse = null; hover = -1; hideTip(); });
    dish.addEventListener('touchmove', e => { const r = dish.getBoundingClientRect(), t = e.touches[0]; mouse = { x: t.clientX - r.left, y: t.clientY - r.top }; }, { passive: true });
    dish.addEventListener('touchend', () => { mouse = null; });
    return {
      start() { size(); if (!running) { running = true; requestAnimationFrame(step); } },
      stop() { running = false; },
      resize: size,
    };
  })();

  /* ================================================================
     MAP: the lab bench of plates
     ================================================================ */
  const bench = $('#bench');
  CH.forEach((c, i) => {
    const b = document.createElement('button');
    b.className = 'plate-btn'; b.dataset.ch = c.id; b.setAttribute('aria-label', 'Plate ' + c.n + ': ' + c.title);
    let cells = '';
    const k = 5 + (i % 3);
    for (let j = 0; j < k; j++) {
      const a = Math.random() * 6.283, d = 6 + Math.random() * 20;
      const ax = (Math.random() - .5) * 26, ay = (Math.random() - .5) * 26, bx = (Math.random() - .5) * 26, by = (Math.random() - .5) * 26;
      cells += `<i style="left:calc(50% + ${Math.cos(a) * d}% - 6px);top:calc(50% + ${Math.sin(a) * d}% - 6px);--s:${10 + Math.random() * 7}px;--c-fill:${c.fill};--c-line:${c.line};--ax:${ax}px;--ay:${ay}px;--bx:${bx}px;--by:${by}px;--dur:${7 + Math.random() * 6}s;--dl:${-Math.random() * 8}s"></i>`;
    }
    b.innerHTML = `<div class="plate-wrap" style="--delay:${-i * 0.8}s"><div class="dish plate"><div class="agar" style="background:radial-gradient(circle at 45% 40%, ${c.fill}, transparent 75%)"></div><span class="num">${c.n}</span>${cells}<div class="glass"></div></div></div>
      <div class="plate-label"><div class="n">PLATE ${c.n}</div><h3>${c.title}</h3><p>${c.teaser}</p></div>`;
    b.addEventListener('click', () => go(c.id, { plate: $('.plate', b) }));
    b.addEventListener('mousemove', e => {
      const r = b.getBoundingClientRect(), x = (e.clientX - r.left) / r.width - .5, y = (e.clientY - r.top) / r.height - .5;
      const p = $('.plate', b); p.style.setProperty('--rx', (-y * 14).toFixed(1) + 'deg'); p.style.setProperty('--ry', (x * 16).toFixed(1) + 'deg');
    });
    b.addEventListener('mouseleave', () => { const p = $('.plate', b); p.style.setProperty('--rx', '0deg'); p.style.setProperty('--ry', '0deg'); });
    bench.appendChild(b);
  });
  const dots = $('#dots');
  CH.forEach((c, i) => {
    const b = document.createElement('button'); b.title = c.n + ' · ' + c.title; b.setAttribute('aria-label', 'Plate ' + c.n + ': ' + c.title);
    b.addEventListener('click', () => go(c.id, { dir: chIdx(cur) >= 0 && i < chIdx(cur) ? -1 : 1 })); dots.appendChild(b);
  });

  /* ================================================================
     NAVIGATION with microscope-style transitions
     ================================================================ */
  let cur = 'home', busy = false;
  const visited = new Set();
  const V = id => document.getElementById('v-' + id);
  const chIdx = id => CH.findIndex(c => c.id === id);
  const DUR = reduce ? 1 : 950, E = 'cubic-bezier(.7,0,.2,1)';
  function anim(el, kf, o) { return el.animate(kf, Object.assign({ duration: DUR, easing: E, fill: 'both' }, o || {})).finished; }
  function viewsRect() { return $('#views').getBoundingClientRect(); }
  function centerOf(el) { const v = viewsRect(), r = el.getBoundingClientRect(); return { x: r.left - v.left + r.width / 2, y: r.top - v.top + r.height / 2, r: r.width / 2, w: r.width }; }
  function maxR() { const v = viewsRect(); return Math.hypot(v.width, v.height); }
  function plateFor(id) { const b = $('.plate-btn[data-ch="' + id + '"]'); return b ? $('.plate', b) : null; }

  async function go(to, opt) {
    opt = opt || {};
    if (busy || to === cur) return;
    busy = true; hideTip(); body.classList.add('hint-off');
    const from = cur, A = V(from), B = V(to);
    const fromCh = chIdx(from) >= 0, toCh = chIdx(to) >= 0;
    B.scrollTop = 0; B.classList.remove('entered');
    B.classList.add('active'); A.style.zIndex = 1; B.style.zIndex = 2;
    beforeEnter(to);
    try {
      if (from === 'home' && to === 'map') await homeToMap(A, B);
      else if (from === 'map' && to === 'home') await mapToHome(A, B);
      else if (from === 'map' && toCh) await lensIn(A, B, centerOf(opt.plate || plateFor(to)), 'zoom');
      else if (fromCh && to === 'map') await lensOut(A, B, plateFor(from));
      else if (fromCh && toCh) {
        const dir = opt.dir || (chIdx(to) > chIdx(from) ? 1 : -1);
        await lensIn(A, B, centerOf(dir > 0 ? $('#next') : $('#prev')), 'slide', dir);
      } else await fade(A, B);
    } catch (e) { /* animation interrupted */ }
    A.classList.remove('active', 'solid', 'entered'); B.classList.remove('solid');
    [A, B].forEach(el => { el.getAnimations().forEach(a => a.cancel()); el.style.zIndex = ''; el.style.transformOrigin = ''; });
    $$('#home-dish, .home-copy, .home-visual .cap').forEach(el => el.getAnimations().forEach(a => a.cancel()));
    cur = to; busy = false;
    if (toCh) visited.add(to);
    chrome(); afterEnter(to);
  }
  async function homeToMap(A, B) {
    const dish = $('#home-dish'), plate = plateFor('question');
    const r1 = dish.getBoundingClientRect(), r2 = plate.getBoundingClientRect();
    const dx = r2.left + r2.width / 2 - (r1.left + r1.width / 2), dy = r2.top + r2.height / 2 - (r1.top + r1.height / 2), s = r2.width / r1.width;
    const c = centerOf(plate); B.style.transformOrigin = c.x + 'px ' + c.y + 'px'; A.style.zIndex = 3;
    await Promise.all([
      anim(dish, [{ transform: 'none' }, { transform: `translate(${dx}px,${dy}px) scale(${s})` }]),
      anim($('.home-copy'), [{ opacity: 1, transform: 'none' }, { opacity: 0, transform: 'translateX(-60px)' }], { duration: DUR * .6 }),
      anim($('.home-visual .cap'), [{ opacity: 1 }, { opacity: 0 }], { duration: DUR * .3 }),
      anim(B, [{ opacity: 0, transform: 'scale(1.35)' }, { opacity: 1, transform: 'scale(1)' }], { delay: DUR * .2, duration: DUR * .8 }),
    ]);
  }
  async function mapToHome(A, B) {
    const dish = $('#home-dish'), plate = plateFor('question');
    const r1 = dish.getBoundingClientRect(), r2 = plate.getBoundingClientRect();
    const dx = r2.left + r2.width / 2 - (r1.left + r1.width / 2), dy = r2.top + r2.height / 2 - (r1.top + r1.height / 2), s = r2.width / r1.width;
    const c = centerOf(plate); A.style.transformOrigin = c.x + 'px ' + c.y + 'px'; B.style.zIndex = 3;
    await Promise.all([
      anim(dish, [{ transform: `translate(${dx}px,${dy}px) scale(${s})` }, { transform: 'none' }]),
      anim($('.home-copy'), [{ opacity: 0, transform: 'translateX(-60px)' }, { opacity: 1, transform: 'none' }], { delay: DUR * .4, duration: DUR * .6 }),
      anim(A, [{ opacity: 1, transform: 'scale(1)' }, { opacity: 0, transform: 'scale(1.35)' }], { duration: DUR * .7 }),
    ]);
  }
  async function lensIn(A, B, c, kind, dir) {
    B.classList.add('solid'); B.style.zIndex = 3;
    const R = maxR(), r0 = kind === 'zoom' ? c.r : 26;
    A.style.transformOrigin = c.x + 'px ' + c.y + 'px';
    await Promise.all([
      anim(B, [{ clipPath: `circle(${r0}px at ${c.x}px ${c.y}px)` }, { clipPath: `circle(${R}px at ${c.x}px ${c.y}px)` }]),
      kind === 'zoom'
        ? anim(A, [{ transform: 'scale(1)', opacity: 1 }, { transform: 'scale(1.7)', opacity: 0 }])
        : anim(A, [{ transform: 'none', opacity: 1 }, { transform: `translateX(${-dir * 70}px) scale(.94)`, opacity: .25 }]),
    ]);
  }
  async function lensOut(A, B, plate) {
    const c = centerOf(plate); A.classList.add('solid'); A.style.zIndex = 3;
    B.style.transformOrigin = c.x + 'px ' + c.y + 'px';
    await Promise.all([
      anim(A, [{ clipPath: `circle(${maxR()}px at ${c.x}px ${c.y}px)` }, { clipPath: `circle(${c.r}px at ${c.x}px ${c.y}px)` }]),
      anim(B, [{ transform: 'scale(1.7)', opacity: 0 }, { transform: 'scale(1)', opacity: 1 }]),
    ]);
  }
  async function fade(A, B) {
    await Promise.all([anim(A, [{ opacity: 1 }, { opacity: 0 }], { duration: DUR * .5 }), anim(B, [{ opacity: 0 }, { opacity: 1 }], { duration: DUR * .6 })]);
  }
  function next() {
    if (cur === 'home') return go('map');
    if (cur === 'map') return go(CH[0].id);
    const i = chIdx(cur); if (i < CH.length - 1) go(CH[i + 1].id, { dir: 1 }); else go('map');
  }
  function prev() {
    if (cur === 'map') return go('home');
    const i = chIdx(cur); if (i === 0) go('map'); else if (i > 0) go(CH[i - 1].id, { dir: -1 });
  }
  function chrome() {
    const i = chIdx(cur);
    body.classList.toggle('in-home', cur === 'home'); body.classList.toggle('in-map', cur === 'map'); body.classList.toggle('in-chapter', i >= 0);
    $$('#dots button').forEach((b, j) => { b.classList.toggle('on', j === i); b.classList.toggle('done', visited.has(CH[j].id) && j !== i); });
    $$('.plate-btn').forEach(b => b.classList.toggle('visited', visited.has(b.dataset.ch)));
    const lab = $('#next .nav-label');
    if (cur === 'home') lab.textContent = 'Explore';
    else if (cur === 'map') lab.textContent = visited.size ? 'Continue' : 'Begin';
    else lab.textContent = i === CH.length - 1 ? 'To the bench' : '';
    $('#next').classList.toggle('labelled', i === CH.length - 1);
  }
  $('#next').addEventListener('click', next);
  $('#prev').addEventListener('click', prev);
  $$('[data-next]').forEach(b => b.addEventListener('click', next));
  $$('[data-go]').forEach(b => b.addEventListener('click', () => go(b.dataset.go)));
  document.addEventListener('keydown', e => {
    if (e.target.tagName === 'SELECT' || e.target.tagName === 'INPUT') return;
    if ($('#drawer').classList.contains('open')) { if (e.key === 'Escape') drawer(false); return; }
    if (e.key === 'ArrowRight') next(); else if (e.key === 'ArrowLeft') prev();
    else if (e.key === 'm' || e.key === 'M' || e.key === 'Escape') { if (cur !== 'map' && cur !== 'home') go('map'); }
  });
  let tx = null, ty = null;
  $('#views').addEventListener('touchstart', e => { tx = e.touches[0].clientX; ty = e.touches[0].clientY; }, { passive: true });
  $('#views').addEventListener('touchend', e => {
    if (tx === null) return; const dx = e.changedTouches[0].clientX - tx, dy = e.changedTouches[0].clientY - ty; tx = null;
    if (e.target.closest('#home-dish')) return;
    if (Math.abs(dx) > 70 && Math.abs(dx) > Math.abs(dy) * 1.5) { dx < 0 ? next() : prev(); }
  }, { passive: true });

  /* ---------------- mode + drawer ---------------- */
  $$('.mode button').forEach(b => b.addEventListener('click', () => {
    $$('.mode button').forEach(x => x.classList.toggle('on', x === b));
    body.classList.toggle('tech', b.dataset.mode === 'tech');
    rerender();
  }));
  function drawer(open) { $('#drawer').classList.toggle('open', open); $('#scrim').classList.toggle('open', open); $('#drawer').setAttribute('aria-hidden', open ? 'false' : 'true'); }
  $('#btn-notes').addEventListener('click', () => drawer(true));
  $('#btn-notes-2').addEventListener('click', () => drawer(true));
  $('#drawer-close').addEventListener('click', () => drawer(false));
  $('#scrim').addEventListener('click', () => drawer(false));
  $$('#drawer-tabs button').forEach(b => b.addEventListener('click', () => {
    $$('#drawer-tabs button').forEach(x => x.classList.toggle('on', x === b));
    $$('.drawer-body .tab').forEach(t => t.classList.toggle('on', t.dataset.t === b.dataset.t));
  }));
  function segGroup(id, cb) {
    const g = $(id); $$('button', g).forEach(b => b.addEventListener('click', () => {
      $$('button', g).forEach(x => x.classList.toggle('on', x === b)); cb(b.dataset.v !== undefined ? b.dataset.v : b.dataset.f);
    }));
  }

  /* ---------------- enter / leave hooks ---------------- */
  let peopleTimer = null;
  function beforeEnter(id) {
    if (id !== 'home') home.stop();
    if (id === 'faith') renderSwarm(true);
    if (id === 'overlap') renderScatter();
    if (id === 'fix') renderVenn();
    if (id === 'external') renderExt();
    if (id === 'question' && !$('#schem svg')) renderSchematic();
  }
  function afterEnter(id) {
    V(id).classList.add('entered');
    if (id === 'home') home.start();
    if (id === 'faith') { countUp($('#faith-pct'), D.faith.median_purity * 100, 0); flySwarm(); }
    if (id === 'overlap') updateScatter(true);
    if (id === 'fix') setTimeout(() => splitVenn(true), 500);
    if (id === 'external') updateExt();
    clearInterval(peopleTimer);
    if (id === 'method') startPeople();
  }
  function rerender() {
    if (cur === 'faith') { renderSwarm(false); }
    if (cur === 'overlap') { renderScatter(); updateScatter(false); }
    if (cur === 'fix') { renderVenn(); splitVenn(false, true); }
    if (cur === 'external') { renderExt(); updateExt(); }
  }
  let rz; window.addEventListener('resize', () => { clearTimeout(rz); rz = setTimeout(() => { home.resize(); rerender(); }, 150); });

  /* ================================================================
     PLATE 01: schematic, built inline
     ================================================================ */
  function renderSchematic() {
    const cat = { a: '#4f7cc0', s: '#8a76c8', g: '#3f8f7a', o: '#c8693a' };
    const cats = ['a', 'a', 'a', 's', 's', 'g', 'g', 'g', 'o', 'o', 'o'];
    const Wd = 700, Hd = 400, ys = cats.map((_, i) => 74 + i * 27), xin = 132, xl = 392, xo = 648;
    const lat = { a: [xl, 140], g: [xl, 278] };
    const plain = !isTech();
    let s = `<svg viewBox="0 0 ${Wd} ${Hd}" role="img" aria-label="Diagram: pathway units can use any input gene but reconstruct only their own pathway genes">`;
    [[xin, 'GENES IN', ''], [(xin + xl) / 2 - 18, 'ENCODER', 'may use any gene'], [xl, plain ? 'PATHWAY UNITS' : 'LATENT UNITS', ''], [(xl + xo) / 2 + 14, 'DECODER', 'own pathway only'], [xo - 10, 'GENES OUT', '']]
      .forEach(([x, a, b]) => { s += `<text x="${x}" y="22" text-anchor="middle" font-size="11.5" letter-spacing="1.4" fill="#5b5f6b">${a}</text><text x="${x}" y="38" text-anchor="middle" font-size="10.5" fill="#8b8a84">${b}</text>`; });
    const encPaths = [];
    cats.forEach((c, i) => {
      const y = ys[i];
      ['a', 'g'].forEach(u => {
        const [lx, ly] = lat[u], leak = c === 'o' || (c === 'g' && u === 'a') || (c === 'a' && u === 'g');
        const d = `M${xin + 8},${y} C${xin + 120},${y} ${lx - 130},${ly} ${lx - 66},${ly}`;
        encPaths.push({ d, c, u });
        s += `<path class="w u-${u} c-${c}${leak ? ' flow-d' : ''}" d="${d}" fill="none" stroke="${cat[c]}" stroke-width="${leak ? 1.3 : 1.5}" opacity="${leak ? .7 : .5}"/>`;
      });
    });
    cats.forEach((c, i) => {
      const y = ys[i];
      ({ a: ['a'], s: ['a', 'g'], g: ['g'], o: [] })[c].forEach(u => {
        const [lx, ly] = lat[u];
        s += `<path class="w u-${u} c-${c}" d="M${lx + 66},${ly} C${lx + 130},${ly} ${xo - 110},${y} ${xo - 8},${y}" fill="none" stroke="${cat[c]}" stroke-width="1.5" opacity=".55"/>`;
      });
    });
    // travelling signals
    if (!reduce) encPaths.filter((_, i) => i % 3 === 0).forEach((p, i) => {
      s += `<circle class="w u-${p.u} c-${p.c}" r="3" fill="${cat[p.c]}"><animateMotion dur="${2.4 + (i % 4) * .5}s" begin="${-i * .37}s" repeatCount="indefinite" path="${p.d}"/></circle>`;
    });
    cats.forEach((c, i) => {
      const y = ys[i];
      s += `<circle cx="${xin}" cy="${y}" r="7.5" fill="#fffdf8" stroke="${cat[c]}" stroke-width="1.6"/><circle cx="${xin + 1.2}" cy="${y - 1}" r="2.8" fill="${cat[c]}" opacity=".75"/>`;
      s += c === 'o' ? `<circle cx="${xo}" cy="${y}" r="7.5" fill="none" stroke="${cat[c]}" stroke-dasharray="2 3" opacity=".5"/>`
        : `<circle cx="${xo}" cy="${y}" r="7.5" fill="#fffdf8" stroke="${cat[c]}" stroke-width="1.6"/><circle cx="${xo + 1.2}" cy="${y - 1}" r="2.8" fill="${cat[c]}" opacity=".75"/>`;
    });
    const labs = { a: ['interferon α', 'genes'], s: ['shared', 'genes'], g: ['interferon γ', 'genes'], o: ['other', 'pathways'] };
    Object.entries({ a: [0, 2], s: [3, 4], g: [5, 7], o: [8, 10] }).forEach(([c, [i0, i1]]) => {
      const y0 = ys[i0] - 8, y1 = ys[i1] + 8, mid = (y0 + y1) / 2;
      s += `<path d="M${xin - 18},${y0} h-5 V${y1} h5" fill="none" stroke="${cat[c]}" stroke-width="1.1"/>`;
      s += `<text x="${xin - 32}" y="${mid - 2}" text-anchor="end" font-size="11.5" fill="${cat[c]}">${labs[c][0]}</text><text x="${xin - 32}" y="${mid + 12}" text-anchor="end" font-size="10.5" fill="${cat[c]}" opacity=".85">${labs[c][1]}</text>`;
    });
    Object.entries(lat).forEach(([u, [lx, ly]]) => {
      const col = cat[u], fill = u === 'a' ? '#e9eff8' : '#e4f0ec';
      s += `<g class="u-node" data-u="${u}"><rect x="${lx - 68}" y="${ly - 25}" width="136" height="50" rx="25" fill="${fill}" stroke="${col}" stroke-width="1.6"/>
        <text x="${lx}" y="${ly - 4}" text-anchor="middle" font-size="10" fill="#8b8a84">${plain ? 'unit named' : 'unit labelled'}</text>
        <text x="${lx}" y="${ly + 12}" text-anchor="middle" font-size="12.5" fill="${col}">${u === 'a' ? 'interferon α' : 'interferon γ'}</text></g>`;
    });
    s += `<text x="${xl}" y="${Hd - 12}" text-anchor="middle" font-size="11" fill="#a9532a">dashed lines: inputs outside the unit’s pathway label</text></svg>`;
    const box = $('#schem'); box.innerHTML = s;
    $$('.u-node', box).forEach(n => {
      n.addEventListener('mouseenter', () => box.setAttribute('data-focus', n.dataset.u));
      n.addEventListener('mouseleave', () => { const on = $('#schem-seg button.on').dataset.f; on ? box.setAttribute('data-focus', on) : box.removeAttribute('data-focus'); });
    });
  }
  segGroup('#schem-seg', f => { const box = $('#schem'); f ? box.setAttribute('data-focus', f) : box.removeAttribute('data-focus'); });

  /* ================================================================
     PLATE 02: rotating held-out person
     ================================================================ */
  $('#people').innerHTML = Array.from({ length: 8 }, () => '<span></span>').join('') ;
  $('#people').insertAdjacentHTML('afterend', '<div class="people-key"><span><i style="background:#3d6fb6"></i>train</span><span><i style="background:#7e6bc4"></i>validation</span><span><i style="background:#c8693a"></i>test (unseen)</span></div>');
  $('#people-new').innerHTML = Array.from({ length: 5 }, () => '<span></span>').join('');
  let pi = 0;
  function startPeople() {
    const ps = $$('#people span');
    const tick = () => { ps.forEach((p, j) => { p.className = j === pi ? 'test' : j === (pi + 1) % 8 ? 'val' : ''; }); pi = (pi + 1) % 8; };
    tick(); if (!reduce) peopleTimer = setInterval(tick, 1300);
  }

  /* ================================================================
     PLATE 03: the swarm of 47 dials + detail card
     ================================================================ */
  let selUnit = least.id, swarmPos = {};
  function renderSwarm(scatterStart) {
    const box = $('#swarm'), W = box.clientWidth || 600, small = W < 560;
    const m = { l: 16, r: 16, t: 34, b: 46 }, r = small ? 7 : 9;
    const x = p => m.l + p * (W - m.l - m.r);
    const placed = []; let cy = 0;
    UNITS.slice().sort((a, b) => a.p - b.p).forEach(u => {
      const px = x(u.p); let k = 0, py = cy;
      for (; k < 60; k++) {
        const off = (k % 2 ? 1 : -1) * Math.ceil(k / 2) * (r * 1.05);
        py = cy + off;
        if (!placed.some(q => Math.hypot(q.x - px, q.y - py) < 2 * r + 1.5)) break;
      }
      placed.push({ x: px, y: py, id: u.id });
    });
    const spread = Math.max(...placed.map(q => Math.abs(q.y))) + r + 8, H = m.t + m.b + 2 * spread;
    placed.forEach(q => { swarmPos[q.id] = [q.x, q.y + m.t + spread]; });
    let s = `<svg viewBox="0 0 ${W} ${H}">`;
    s += `<rect x="${x(0)}" y="${m.t - 6}" width="${x(.5) - x(0)}" height="${H - m.t - m.b + 12}" rx="14" fill="#f8e5d6" opacity=".55"/>`;
    s += `<text class="anno-s" x="${x(.25)}" y="${m.t - 14}" text-anchor="middle">${isTech() ? 'purity < 0.5' : 'mostly outside its pathway'}</text><text class="anno-s" x="${x(.75)}" y="${m.t - 14}" text-anchor="middle">${isTech() ? 'purity ≥ 0.5' : 'at least half on its own pathway'}</text>`;
    s += `<line class="ax" x1="${m.l}" x2="${W - m.r}" y1="${H - m.b + 10}" y2="${H - m.b + 10}"/>`;
    [0, .25, .5, .75, 1].forEach(v => { s += `<text class="tick" x="${x(v)}" y="${H - m.b + 26}" text-anchor="middle">${isTech() ? v.toFixed(2) : Math.round(v * 100) + '%'}</text>`; });
    s += `<text class="axlab" x="${(m.l + W - m.r) / 2}" y="${H - 4}" text-anchor="middle">${isTech() ? 'median attribution purity (absolute IG share on member genes)' : 'share of attribution on its own pathway genes'}</text>`;
    s += `<line x1="${x(.5)}" x2="${x(.5)}" y1="${m.t - 6}" y2="${H - m.b + 10}" stroke="#c8693a" stroke-dasharray="4 4"/>`;
    s += `<circle class="ring" r="${r + 5}" style="transform:translate(${swarmPos[selUnit][0]}px,${swarmPos[selUnit][1]}px)"/>`;
    UNITS.forEach((u, i) => {
      const [px, py] = swarmPos[u.id], low = u.p < .5;
      const start = scatterStart && !reduce ? `translate(${Math.random() * W}px,${-40 - Math.random() * 60}px)` : `translate(${px}px,${py}px)`;
      s += `<g class="cellg mv" data-id="${u.id}" style="transform:${start};transition-delay:${(i * 14)}ms">
        <circle class="mem" r="${r}" fill="${low ? '#f3cdb3' : '#c4d7ef'}" stroke="${low ? '#c8693a' : '#3d6fb6'}" stroke-width="1.4"/>
        <circle r="${r * .36}" cx="${r * .2}" cy="${-r * .15}" fill="${low ? '#c8693a' : '#3d6fb6'}" opacity=".55" pointer-events="none"/></g>`;
    });
    s += '</svg>'; box.innerHTML = s;
    $$('.cellg', box).forEach(g => {
      const u = byId[g.dataset.id];
      g.addEventListener('mouseenter', e => showTip(`<b>${u.n}</b>${isTech() ? 'purity ' + f3(u.p) + ' · Spearman ' + f3(u.s) + ' · ablation ' + u.a.toFixed(1) + '×' : pct(u.p) + ' of attribution on its own pathway genes'}`, e));
      g.addEventListener('mousemove', e => showTip(tip.innerHTML, e));
      g.addEventListener('mouseleave', hideTip);
      g.addEventListener('click', () => selectUnit(u.id));
    });
    renderUnitCard();
  }
  function flySwarm() {
    requestAnimationFrame(() => requestAnimationFrame(() => {
      $$('#swarm .cellg').forEach(g => { const [px, py] = swarmPos[g.dataset.id]; g.style.transform = `translate(${px}px,${py}px)`; });
    }));
  }
  function selectUnit(id) {
    selUnit = id; const [px, py] = swarmPos[id];
    const ring = $('#swarm .ring'); if (ring) ring.style.transform = `translate(${px}px,${py}px)`;
    renderUnitCard();
  }
  function renderUnitCard() {
    const u = byId[selUnit], nOut = u.top.filter(g => !g[2]).length, mx = Math.max(...u.top.map(g => g[1]));
    const opts = UNITS.slice().sort((a, b) => a.n.localeCompare(b.n)).map(x => `<option value="${x.id}" ${x.id === u.id ? 'selected' : ''}>${x.n}</option>`).join('');
    const say = nOut === 0 ? 'All 10 of its top 10 attributed genes are in its named pathway.'
      : `${nOut} of its top 10 attributed genes are <b>not</b> in its named pathway.`;
    $('#unit-card').innerHTML = `<div>
        <div class="who"><small>${u.id === least.id ? 'lowest-purity unit' : u.id === purest.id ? 'highest-purity unit' : 'pathway unit'}</small>${u.n}</div>
        <div class="share">${pct(u.p)}<small>${isTech() ? 'attribution purity ' + f3(u.p) : 'of its attribution on its own pathway genes'}</small></div>
        <p class="say p">${say}</p>
        <p class="say t">Spearman ${f3(u.s)} (min fold ${f3(u.smin)}) · ablation ratio ${u.a.toFixed(1)}× (min ${u.amin.toFixed(1)}) · enrichment ${u.en.toFixed(1)}× · ${u.size} genes</p>
        <select id="unit-sel" aria-label="Choose a pathway unit">${opts}</select>
        <div class="seg" style="margin-top:.6rem"><button data-u="${least.id}">Lowest purity</button><button data-u="${IFNA}">Interferon α</button><button data-u="${purest.id}">Highest purity</button></div>
      </div>
      <div><div class="chart-title">Top 10 input genes by mean |IG|</div><div class="genebars">${u.top.map(g => `<div class="gb" title="${g[0]}: mean |IG| ${g[1]}"><span>${g[0]}</span><span class="tr"><i class="${g[2] ? 'm' : 'o'}" data-w="${(g[1] / mx * 100).toFixed(1)}"></i></span></div>`).join('')}</div>
      <div class="legend"><span><i style="background:rgba(61,111,182,.82)"></i>in its named pathway</span><span><i style="background:#d9cfbb"></i>outside its pathway</span></div></div>`;
    $('#unit-sel').addEventListener('change', e => selectUnit(e.target.value));
    $$('#unit-card .seg button').forEach(b => { b.classList.toggle('on', b.dataset.u === u.id); b.addEventListener('click', () => selectUnit(b.dataset.u)); });
    requestAnimationFrame(() => requestAnimationFrame(() => $$('#unit-card .tr i').forEach(i => { i.style.width = i.dataset.w + '%'; })));
  }

  /* ================================================================
     PLATE 04: overlap scatter (dots glide between states)
     ================================================================ */
  const OV = D.overlap; let ovCells = 'all', ovModel = 'D', ovGeom = null;
  function renderScatter() {
    const box = $('#scatter'), W = box.clientWidth || 640, H = clamp(Math.round(W * .66), 290, 450), m = { l: 48, r: 16, t: 16, b: 48 };
    const jmax = Math.max(...OV.pairs.map(p => p[2])) * 1.06;
    const x = j => m.l + j / jmax * (W - m.l - m.r), y = c => H - m.b - c * (H - m.t - m.b);
    ovGeom = { x, y, W, H, m };
    let s = `<svg viewBox="0 0 ${W} ${H}">`;
    [0, .25, .5, .75, 1].forEach(v => { s += `<line class="grid-l" x1="${m.l}" x2="${W - m.r}" y1="${y(v)}" y2="${y(v)}"/><text class="tick" x="${m.l - 8}" y="${y(v) + 3}" text-anchor="end">${v}</text>`; });
    [0, .1, .2, .3].forEach(v => { s += `<text class="tick" x="${x(v)}" y="${H - m.b + 16}" text-anchor="middle">${isTech() ? v.toFixed(1) : Math.round(v * 100) + '%'}</text>`; });
    s += `<text class="axlab" x="${(m.l + W - m.r) / 2}" y="${H - 8}" text-anchor="middle">${isTech() ? 'Jaccard overlap of the two fold-level gene sets' : 'gene-list overlap (shared genes ÷ total unique genes)'}</text>`;
    s += `<text class="axlab" transform="translate(13 ${(m.t + H - m.b) / 2}) rotate(-90)" text-anchor="middle">${isTech() ? '|Spearman| between unit scores' : 'how similar the two units’ scores are'}</text>`;
    s += `<path id="ov-trend" fill="none" stroke="#23262f" stroke-width="2" stroke-linecap="round" opacity=".7" style="transition:d 1s cubic-bezier(.2,.8,.2,1)"/>`;
    s += `<g id="ov-dots">`;
    OV.pairs.forEach((p, i) => {
      if (i === OV.ifn) return;
      s += `<circle class="mv" r="3.1" data-i="${i}" fill="rgba(61,111,182,.42)" style="transform:translate(${x(p[2])}px,${y(0)}px);transition-delay:${(i % 40) * 6}ms"/>`;
    });
    s += `</g><g id="ov-ifn-g" class="mv" style="transform:translate(${x(OV.pairs[OV.ifn][2])}px,${y(0)}px)"><circle r="8.5" fill="#c8693a" stroke="#fffdf8" stroke-width="2"/>
      <text class="anno" x="-14" y="-14" text-anchor="end">interferon α × γ</text></g>`;
    s += `<text id="ov-trend-lab" class="anno-s" x="${W - m.r}" y="${m.t + 12}" text-anchor="end">— median correlation within overlap bands</text></svg>`;
    box.innerHTML = s;
    $$('#ov-dots circle', box).forEach(c => {
      c.addEventListener('mouseenter', e => { const p = OV.pairs[+c.dataset.i], v = OV.corr[ovModel + '|' + ovCells][+c.dataset.i]; showTip(`<b>${OV.names[p[0]]} × ${OV.names[p[1]]}</b>${isTech() ? 'Jaccard ' + p[2].toFixed(3) + ' · |ρ| ' + f3(v) : Math.round(p[2] * 100) + '% gene-set overlap · correlation ' + f3(v)}`, e); });
      c.addEventListener('mouseleave', hideTip);
    });
  }
  function trendPath(vals) {
    const edges = [[0, 0], [1e-9, .02], [.02, .04], [.04, .07], [.07, .12], [.12, .35]];
    const pts = [];
    edges.forEach(([lo, hi], k) => {
      const js = [], cs = [];
      OV.pairs.forEach((p, i) => { const j = p[2]; if ((k === 0 ? j === 0 : j > lo && j <= hi)) { js.push(j); cs.push(vals[i]); } });
      if (cs.length < 8) return;
      const med = a => { const b = a.slice().sort((u, v) => u - v), n = b.length; return n % 2 ? b[(n - 1) / 2] : (b[n / 2 - 1] + b[n / 2]) / 2; };
      pts.push([ovGeom.x(med(js)), ovGeom.y(med(cs))]);
    });
    return 'M' + pts.map(p => p[0].toFixed(1) + ',' + p[1].toFixed(1)).join(' L');
  }
  function updateScatter(animate) {
    if (!ovGeom) return;
    const key = ovModel + '|' + ovCells, vals = OV.corr[key], S = OV.summary[key], { x, y } = ovGeom;
    $$('#ov-dots circle').forEach(c => { const i = +c.dataset.i; c.style.transform = `translate(${x(OV.pairs[i][2])}px,${y(vals[i])}px)`; });
    $('#ov-ifn-g').style.transform = `translate(${x(OV.pairs[OV.ifn][2])}px,${y(vals[OV.ifn])}px)`;
    $('#ov-trend').setAttribute('d', trendPath(vals));
    const el = $('#ov-ifn'); animate ? countUp(el, S.ifn, 2, 900) : (el.textContent = S.ifn.toFixed(2), el.dataset.cur = S.ifn);
    const Sa = OV.summary[ovModel + '|all'], Sc = OV.summary[ovModel + '|ctrl_only'];
    $('#ov-ifn-sub').textContent = isTech() ? `|ρ| = ${f3(Sa.ifn)} across all cells; |ρ| = ${f3(Sc.ifn)} in CTRL-only cells` : (ovCells === 'all' ? 'strong co-movement across all cells' : 'little correlation in resting cells');
    $('#ov-stats').innerHTML = `<span>ρ(overlap, |corr|)</span><span>${f3(S.rho)}</span><span>median |corr|, pairs that share genes (${S.n_ov})</span><span>${f3(S.ov)}</span><span>median |corr|, pairs that share none (${S.n - S.n_ov})</span><span>${f3(S.non)}</span>`;
  }
  segGroup('#ov-cells', v => { ovCells = v; updateScatter(true); });
  segGroup('#ov-model', v => { ovModel = v; updateScatter(true); });

  /* ================================================================
     PLATE 05: the Venn that splits apart
     ================================================================ */
  const NA = D.ifn.IFNA_UNIQUE.length, NS = D.ifn.IFN_SHARED.length, NG = D.ifn.IFNG_UNIQUE.length;
  let venn = null, vennState = 0;
  function renderVenn() {
    const box = $('#venn'), W = box.clientWidth || 600, H = clamp(Math.round(W * .56), 260, 380);
    const sc = Math.min(W, H * 1.6) / 62, ra = sc * Math.sqrt(NA + NS), rg = sc * Math.sqrt(NG + NS);
    const target = Math.PI * sc * sc * NS;
    const lens = d => {
      if (d >= ra + rg) return 0; if (d <= Math.abs(ra - rg)) return Math.PI * Math.min(ra, rg) ** 2;
      const a1 = ra * ra * Math.acos((d * d + ra * ra - rg * rg) / (2 * d * ra)), a2 = rg * rg * Math.acos((d * d + rg * rg - ra * ra) / (2 * d * rg));
      return a1 + a2 - .5 * Math.sqrt((-d + ra + rg) * (d + ra - rg) * (d - ra + rg) * (d + ra + rg));
    };
    let lo = Math.abs(ra - rg), hi = ra + rg; for (let k = 0; k < 60; k++) { const mid = (lo + hi) / 2; if (lens(mid) > target) lo = mid; else hi = mid; }
    const d = (lo + hi) / 2, cy = H * .5, cxA = W / 2 - (d + rg - ra) / 2 - ra * .1, cxG = cxA + d;
    const r1 = sc * Math.sqrt(NA), r2 = sc * Math.sqrt(NS), r3 = sc * Math.sqrt(NG), gap = W * .04;
    const tot = 2 * (r1 + r2 + r3) + 2 * gap, x1 = (W - tot) / 2 + r1, x2 = x1 + r1 + gap + r2, x3 = x2 + r2 + gap + r3;
    const lensX = (cxA + ra + cxG - rg) / 2;
    venn = {
      A: { a: [cxA, cy, ra], b: [x1, cy, r1] }, S: { a: [lensX, cy, 0.1], b: [x2, cy, r2] }, G: { a: [cxG, cy, rg], b: [x3, cy, r3] },
      tA: { a: [cxA - ra * .5, cy], b: [x1, cy] }, tS: { a: [lensX, cy], b: [x2, cy] }, tG: { a: [cxG + rg * .4, cy], b: [x3, cy] }, W, H,
    };
    const T = isTech();
    box.innerHTML = `<svg viewBox="0 0 ${W} ${H}"><defs><clipPath id="clipA"><circle id="clipAc" cx="${cxA}" cy="${cy}" r="${ra}"/></clipPath></defs>
      <circle id="vG" cx="${cxG}" cy="${cy}" r="${rg}" fill="#e4f0ec" stroke="#3f8f7a" stroke-width="1.6"/>
      <circle id="vA" cx="${cxA}" cy="${cy}" r="${ra}" fill="#dfe8f6" fill-opacity=".85" stroke="#4f7cc0" stroke-width="1.6"/>
      <circle id="vLens" cx="${cxG}" cy="${cy}" r="${rg}" fill="#d6cdf0" clip-path="url(#clipA)" style="transition:opacity .5s"/>
      <circle id="vS" cx="${lensX}" cy="${cy}" r="0.1" fill="#d6cdf0" stroke="#8a76c8" stroke-width="1.6"/>
      <text id="tA" class="vnum" x="${venn.tA.a[0]}" y="${cy + 9}" text-anchor="middle" font-family="Instrument Serif, Georgia, serif" font-size="28" fill="#23262f">${NA}</text>
      <text id="tS" class="vnum" x="${lensX}" y="${cy + 11}" text-anchor="middle" font-family="Instrument Serif, Georgia, serif" font-size="34" fill="#23262f">${NS}</text>
      <text id="tG" class="vnum" x="${venn.tG.a[0]}" y="${cy + 9}" text-anchor="middle" font-family="Instrument Serif, Georgia, serif" font-size="28" fill="#23262f">${NG}</text>
      <g id="labA" class="vlabel"><text x="${cxA - ra * .2}" y="${cy - ra - 14}" text-anchor="middle" class="anno" fill="#4f7cc0">${T ? 'Hallmark IFN-α · ' + (NA + NS) : 'interferon α list · ' + (NA + NS) + ' genes'}</text></g>
      <g id="labG" class="vlabel"><text x="${cxG + rg * .2}" y="${cy + rg + 22}" text-anchor="middle" class="anno" fill="#3f8f7a">${T ? 'Hallmark IFN-γ · ' + (NG + NS) : 'interferon γ list · ' + (NG + NS) + ' genes'}</text></g>
      <g id="labSplit" class="vlabel" style="opacity:0">
        <text x="${x1}" y="${cy + Math.max(r1, r3) + 26}" text-anchor="middle" class="anno-s">${T ? 'IFNA_UNIQUE' : 'α only'}</text>
        <text x="${x2}" y="${cy + Math.max(r1, r3) + 26}" text-anchor="middle" class="anno">${T ? 'IFN_SHARED' : 'shared → own unit'}</text>
        <text x="${x3}" y="${cy + Math.max(r1, r3) + 26}" text-anchor="middle" class="anno-s">${T ? 'IFNG_UNIQUE' : 'γ only'}</text></g></svg>`;
    vennState = 0;
  }
  function tweenVenn(k) {
    const L = (p, q) => p + (q - p) * k, set = (id, a) => { const e = document.getElementById(id); e.setAttribute('cx', a[0]); e.setAttribute('cy', a[1]); e.setAttribute('r', Math.max(.1, a[2])); };
    const P = o => [L(o.a[0], o.b[0]), L(o.a[1], o.b[1]), L(o.a[2], o.b[2])];
    set('vA', P(venn.A)); set('clipAc', P(venn.A)); set('vS', P(venn.S)); set('vG', P(venn.G)); set('vLens', P(venn.G));
    ['tA', 'tS', 'tG'].forEach(id => { const o = venn[id], e = document.getElementById(id); e.setAttribute('x', L(o.a[0], o.b[0])); });
    $('#vLens').style.opacity = String(1 - k);
    $('#labA').style.opacity = $('#labG').style.opacity = String(1 - k);
    $('#labSplit').style.opacity = String(k);
  }
  function splitVenn(animate, instant) {
    if (!venn) return;
    const to = 1, from = vennState;
    if (instant || reduce) { tweenVenn(to); vennState = to; showGenes(); return; }
    const t0 = performance.now(), dur = 1500;
    (function step(t) { const k = Math.min(1, (t - t0) / dur); tweenVenn(from + (to - from) * easeIO(k)); if (k < 1) requestAnimationFrame(step); else { vennState = to; showGenes(); } })(t0);
  }
  function showGenes() {
    const g = D.ifn.IFN_SHARED, n = 28;
    $('#venn-genes').innerHTML = `<span style="border:0;background:none;padding-left:0;opacity:1;animation:none">${isTech() ? 'IFN_SHARED genes:' : 'The 72 shared genes:'}</span>` +
      g.slice(0, n).map((x, i) => `<span style="animation-delay:${i * 25}ms">${x}</span>`).join('') + `<span style="animation-delay:${n * 25}ms">+${g.length - n} more</span>`;
  }
  $('#fix-replay').addEventListener('click', () => { if (!venn) return; tweenVenn(0); vennState = 0; $('#venn-genes').innerHTML = ''; setTimeout(() => splitVenn(true), 350); });

  /* ================================================================
     PLATE 06: external validation
     ================================================================ */
  let extCell = 'CD14 Mono', extStim = 'IFNa', extGeom = null;
  const REPS = ['A', 'D', 'D2', 'A2'];
  const REP_PLAIN = { A: ['Simple average', 'original sets'], D: ['Neural model', 'original sets'], D2: ['Neural model', 'split sets'], A2: ['Simple average', 'split sets'] };
  const REP_TECH = { A: ['A', 'set average'], D: ['D', 'pathway AE'], D2: ['D2', 'overlap-aware AE'], A2: ['A2', 'split average'] };
  const REP_COL = { A: ['#dfe8f6', '#3d6fb6'], D: ['#3d6fb6', '#3d6fb6'], D2: ['#c8693a', '#c8693a'], A2: ['#f8e5d6', '#c8693a'] };
  function renderExt() {
    const box = $('#ext'), W = box.clientWidth || 640, H = clamp(Math.round(W * .62), 300, 430), m = { l: 44, r: 12, t: 24, b: 58 };
    const rows = D.sep.filter(r => r[0] === extStim), ymax = Math.ceil(Math.max(...rows.map(r => r[4])) + .3);
    const band = (W - m.l - m.r) / 4, bw = Math.min(70, band * .52), y = v => H - m.b - v / ymax * (H - m.t - m.b);
    extGeom = { W, H, m, band, bw, y, ymax };
    const T = isTech(), lab = T ? REP_TECH : REP_PLAIN;
    let s = `<svg viewBox="0 0 ${W} ${H}">`;
    for (let v = 0; v <= ymax; v++) s += `<line class="grid-l" x1="${m.l}" x2="${W - m.r}" y1="${y(v)}" y2="${y(v)}"/><text class="tick" x="${m.l - 8}" y="${y(v) + 3}" text-anchor="end">${v}</text>`;
    s += `<text class="axlab" transform="translate(12 ${(m.t + H - m.b) / 2}) rotate(-90)" text-anchor="middle">${T ? 'SEP (Kang-SD units)' : 'separation (higher = clearer)'}</text>`;
    REPS.forEach((r, i) => {
      const cx = m.l + band * (i + .5), [fill, stroke] = REP_COL[r];
      s += `<rect class="bar mv" data-r="${r}" x="${cx - bw / 2}" y="${m.t}" width="${bw}" height="${H - m.b - m.t}" rx="8" fill="${fill}" stroke="${stroke}" stroke-width="${r === 'A' || r === 'A2' ? 2 : 0}"
        style="transform-box:fill-box;transform-origin:50% 100%;transform:scaleY(0.001)"/>`;
      s += `<text class="anno med" data-r="${r}" x="${cx + bw / 2 + 6}" y="${H - m.b - 4}" text-anchor="start" style="transition:transform 1s cubic-bezier(.2,.8,.2,1);font-family:var(--mono);font-size:11.5px">–</text>`;
      s += `<text class="anno" x="${cx}" y="${H - m.b + 20}" text-anchor="middle">${lab[r][0]}</text><text class="anno-s" x="${cx}" y="${H - m.b + 36}" text-anchor="middle">${lab[r][1]}</text>`;
    });
    s += `<g id="ext-dots"></g></svg>`;
    box.innerHTML = s;
    renderScore();
  }
  function updateExt() {
    if (!extGeom) return;
    const { H, m, band, y } = extGeom, rows = D.sep.filter(r => r[0] === extStim && r[1] === extCell);
    const donors = [...new Set(rows.map(r => r[3]))].sort();
    const med = a => { const b = a.slice().sort((p, q) => p - q), n = b.length; return n % 2 ? b[(n - 1) / 2] : (b[n / 2 - 1] + b[n / 2]) / 2; };
    let best = null, bestv = -1e9;
    REPS.forEach(r => {
      const vals = rows.filter(x => x[2] === r).map(x => x[4]), md = med(vals);
      if (md > bestv) { bestv = md; best = r; }
      const bar = $(`#ext rect.bar[data-r="${r}"]`); bar.style.transform = `scaleY(${Math.max(.002, md / extGeom.ymax)})`;
      const t = $(`#ext text.med[data-r="${r}"]`); t.textContent = (isTech() ? 'median ' : '') + md.toFixed(2); t.style.transform = `translateY(${y(md) - (H - m.b) + 4}px)`;
    });
    const g = $('#ext-dots'), offs = [-14, -7, 0, 7, 14];
    if (!g.children.length) {
      let s = '';
      REPS.forEach((r, i) => donors.forEach((dn, k) => {
        s += `<circle class="mv dd" data-r="${r}" data-k="${k}" r="4.2" fill="#fffdf8" stroke="#23262f" stroke-width="1.3" style="transform:translate(${m.l + band * (i + .5) + offs[k]}px,${H - m.b}px)"/>`;
      }));
      g.innerHTML = s;
      $$('circle.dd', g).forEach(c => {
        c.addEventListener('mouseenter', e => { const row = D.sep.find(x => x[0] === extStim && x[1] === extCell && x[2] === c.dataset.r && x[3] === donors[+c.dataset.k]); if (row) showTip(`<b>Donor ${row[3]}</b>${(isTech() ? REP_TECH : REP_PLAIN)[row[2]].join(' · ')}: ${row[4].toFixed(2)}`, e); });
        c.addEventListener('mouseleave', hideTip);
      });
    }
    requestAnimationFrame(() => $$('circle.dd', g).forEach(c => {
      const r = c.dataset.r, k = +c.dataset.k, row = rows.find(x => x[2] === r && x[3] === donors[k]), i = REPS.indexOf(r);
      if (row) c.style.transform = `translate(${m.l + band * (i + .5) + offs[k]}px,${y(row[4])}px)`;
    }));
    const cellName = { 'CD14 Mono': 'monocytes', B: 'B cells', NK: 'NK cells', T: 'T cells' }[extCell];
    $('#ext-title').textContent = isTech() ? `SEP per donor · ${extCell} · ${extStim === 'IFNa' ? 'IFN-α' : 'IFN-β'} vs IFN-γ` : `Type-I vs type-II interferon separation · higher is better · ${cellName}`;
    renderScore();
  }
  function renderScore() {
    const T = isTech(), cells = ['CD14 Mono', 'B', 'NK', 'T'], short = { 'CD14 Mono': 'Mono', B: 'B', NK: 'NK', T: 'T' };
    const rows = [
      ['A2>A', T ? 'A2 vs A · decomposition-only comparison' : 'Split average beats original average'],
      ['D2>D', T ? 'D2 vs D · primary H5 (CD14, IFN-α)' : 'Split neural model beats original neural model'],
      ['D2>A2', T ? 'D2 vs A2 · does the neural model add beyond the split?' : 'Split neural model beats split average'],
    ];
    if (T) rows.push(['D2>A', 'D2 vs A']);
    let s = `<div class="sc-head"><div>${T ? 'Comparison (≥4/5 donors & median > 0)' : 'On each cell type'}</div><div>${cells.map(c => `<span>${short[c]}</span>`).join('')}</div></div>`;
    rows.forEach(([k, q]) => {
      s += `<div class="sc-row"><div class="q">${q}</div><div class="sc-cells">` + cells.map(c => {
        const v = D.verdicts[`${extStim}|${c}|${k}`];
        return `<span class="${v.met ? 'met' : 'no'}${c === extCell ? ' cur' : ''}" title="${v.nb}/${v.n} donors · median ${sgn(v.med)}"><b>${v.met ? '✓' : '✕'}</b>${v.nb}/${v.n}</span>`;
      }).join('') + `</div></div>`;
    });
    if (T) {
      const c = 'CD14 Mono', a2 = D.verdicts[`${extStim}|${c}|A2>A`], d2 = D.verdicts[`${extStim}|${c}|D2>D`], cr = D.verdicts[`${extStim}|${c}|split>net`];
      s += `<div class="small">CD14 credit rule: A2−A = ${sgn(a2.med)} vs D2−D = ${sgn(d2.med)}; Δsplit−Δnet = ${sgn(cr.med)} → improvement credited primarily to the gene-set decomposition.`
        + (extCell !== c ? ` ${extCell} is a secondary replication and does not redefine the primary H5 verdict.` : '') + `</div>`;
    }
    $('#scorecard').innerHTML = s;
  }
  segGroup('#ext-cell', v => { extCell = v; updateExt(); });
  segGroup('#ext-stim', v => { extStim = v; renderExt(); updateExt(); });

  /* ---------------- start ---------------- */
  $$('[data-count]').forEach(e => countUp(e, +e.dataset.count, 0, 1400));
  chrome(); home.start();
  setTimeout(() => body.classList.add('hint-off'), 9000);
})();
