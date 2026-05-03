// Copy-to-clipboard — works for both .code-block and .deriv-code-panel
document.querySelectorAll('.copy-btn').forEach(btn => {
  btn.addEventListener('click', () => {
    let text;
    const block = btn.closest('.code-block');
    const panel = btn.closest('.deriv-code-panel');

    if (block) {
      text = block.querySelector('pre code').innerText;
    } else if (panel) {
      text = panel.querySelector('pre').innerText;
    }
    if (!text) return;

    navigator.clipboard.writeText(text).then(() => {
      btn.textContent = 'Copied!';
      btn.classList.add('copied');
      setTimeout(() => { btn.textContent = 'Copy'; btn.classList.remove('copied'); }, 2000);
    });
  });
});

// Highlight active sidebar link on scroll
const sections = document.querySelectorAll('[id]');
const navLinks = document.querySelectorAll('#sidebar nav a[href^="#"]');

const observer = new IntersectionObserver(entries => {
  entries.forEach(entry => {
    if (entry.isIntersecting) {
      navLinks.forEach(a => a.classList.remove('active'));
      const active = document.querySelector(`#sidebar nav a[href="#${entry.target.id}"]`);
      if (active) active.classList.add('active');
    }
  });
}, { rootMargin: '-20% 0px -70% 0px' });

sections.forEach(s => observer.observe(s));
