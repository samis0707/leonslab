export const languages = {
  de: 'DE',
  en: 'EN',
} as const;

export const defaultLang = 'de';
export type Lang = keyof typeof languages;

export const ui = {
  de: {
    'nav.cartoons': 'Cartoons',
    'nav.stickers': 'Sticker',
    'nav.labSolutions': 'Lab Solutions',
    'nav.park': 'Pseudomonas Park',
    'nav.ideas': 'Ideen',

    'hero.eyebrow': 'LeonsLab',
    'hero.title.line1': 'Wo Mikrobiologie',
    'hero.title.line2': 'auf Gestaltung trifft.',
    'hero.subtitle': 'Ein Creative Space von Leon Rösch — Forschung, Cartoons, Sticker und Werkstatt-Lösungen aus dem Labor.',
    'hero.cta.work': 'Arbeiten ansehen',
    'hero.cta.park': 'Pseudomonas Park betreten',

    'about.eyebrow': 'Über',
    'about.title': 'Naturwissenschaft mit Handschrift.',
    'about.body': 'Ich bin Leon — Mikrobiologe, Zeichner und Bastler. Hier sammle ich, was zwischen Petrischale und Skizzenblock entsteht: ernsthafte Forschung neben verspielten Cartoons, durchdachten Sticker-Designs und improvisierten Lösungen aus dem Laboralltag.',

    'sections.work': 'Arbeit',
    'sections.exploring': 'In Arbeit',
    'sections.viewAll': 'Alle ansehen',

    'cards.cartoons.title': 'Cartoons',
    'cards.cartoons.desc': 'Mikrobiologie mit Augenzwinkern — eine Serie zusammenhängender Zeichnungen.',
    'cards.stickers.title': 'Sticker',
    'cards.stickers.desc': 'Designs, die ihren Weg auf Notebooks, Wände und Beton finden.',
    'cards.labSolutions.title': 'Lab Solutions',
    'cards.labSolutions.desc': 'Improvisierte Lösungen für reale Probleme im Labor — designed & gefertigt.',

    'park.eyebrow': 'Creative Space',
    'park.title': 'Pseudomonas Park',
    'park.tagline': 'Der Ort, an dem ich als Künstler arbeite.',
    'park.body': 'Pseudomonas Park ist mehr als ein Sticker. Es ist die Welt, in der meine kreativen Projekte zuhause sind — ein Studio, ein Labor, ein Spielfeld.',

    'footer.contact': 'Kontakt',
    'footer.colophon': 'Gebaut mit Astro · System-Schrift',

    'lang.switch': 'EN',
    'lang.switch.aria': 'Sprache wechseln zu Englisch',

    'page.cartoons.title': 'Cartoons',
    'page.cartoons.subtitle': 'Eine wachsende Sammlung mikrobiologischer Zeichnungen.',
    'page.cartoons.empty': 'Bilder folgen in Kürze.',
    'gallery.filter.all': 'Alle',
    'gallery.tag.lab-life': 'Lab Life',
    'gallery.tag.slice-of-life': 'Slice of Life',
    'gallery.tag.pop-parody': 'Pop-Parodie',
    'gallery.tag.pseudomonas-park': 'Pseudomonas Park',
    'gallery.count': 'Cartoons',
    'gallery.empty': 'Keine Cartoons in dieser Kategorie.',
    'lightbox.close': 'Schließen',
    'lightbox.prev': 'Vorheriges Bild',
    'lightbox.next': 'Nächstes Bild',
    'page.stickers.title': 'Sticker',
    'page.stickers.subtitle': 'Designs aus dem Studio.',
    'page.labSolutions.title': 'Lab Solutions',
    'page.labSolutions.subtitle': 'Wenn ein Teil fehlt, baut man es.',
    'page.park.title': 'Pseudomonas Park',
    'page.park.subtitle': 'Der Creative Space.',
    'page.ideas.title': 'Ideen',
    'page.ideas.subtitle': 'Halbfertiges, Wildes, Brauendes.',
  },
  en: {
    'nav.cartoons': 'Cartoons',
    'nav.stickers': 'Stickers',
    'nav.labSolutions': 'Lab Solutions',
    'nav.park': 'Pseudomonas Park',
    'nav.ideas': 'Ideas',

    'hero.eyebrow': 'LeonsLab',
    'hero.title.line1': 'Where microbiology',
    'hero.title.line2': 'meets making.',
    'hero.subtitle': 'A creative space by Leon Rösch — research, cartoons, stickers and workshop solutions from the lab.',
    'hero.cta.work': 'See the work',
    'hero.cta.park': 'Enter Pseudomonas Park',

    'about.eyebrow': 'About',
    'about.title': 'Science with a personal hand.',
    'about.body': "I'm Leon — microbiologist, illustrator and tinkerer. This is where I collect what happens between petri dish and sketchbook: serious research alongside playful cartoons, considered sticker designs and improvised solutions for everyday lab problems.",

    'sections.work': 'Work',
    'sections.exploring': 'In progress',
    'sections.viewAll': 'View all',

    'cards.cartoons.title': 'Cartoons',
    'cards.cartoons.desc': 'Microbiology with a wink — a series of connected drawings.',
    'cards.stickers.title': 'Stickers',
    'cards.stickers.desc': 'Designs that find their way onto laptops, walls and concrete.',
    'cards.labSolutions.title': 'Lab Solutions',
    'cards.labSolutions.desc': 'Improvised answers to real problems in the lab — designed & fabricated.',

    'park.eyebrow': 'Creative Space',
    'park.title': 'Pseudomonas Park',
    'park.tagline': 'The place where I work as an artist.',
    'park.body': "Pseudomonas Park is more than a sticker. It's the world my creative projects live in — a studio, a lab, a playground.",

    'footer.contact': 'Contact',
    'footer.colophon': 'Built with Astro · System type',

    'lang.switch': 'DE',
    'lang.switch.aria': 'Switch language to German',

    'page.cartoons.title': 'Cartoons',
    'page.cartoons.subtitle': 'A growing collection of microbiology drawings.',
    'page.cartoons.empty': 'Images coming soon.',
    'gallery.filter.all': 'All',
    'gallery.tag.lab-life': 'Lab Life',
    'gallery.tag.slice-of-life': 'Slice of Life',
    'gallery.tag.pop-parody': 'Pop Parody',
    'gallery.tag.pseudomonas-park': 'Pseudomonas Park',
    'gallery.count': 'cartoons',
    'gallery.empty': 'No cartoons in this category.',
    'lightbox.close': 'Close',
    'lightbox.prev': 'Previous image',
    'lightbox.next': 'Next image',
    'page.stickers.title': 'Stickers',
    'page.stickers.subtitle': 'Designs from the studio.',
    'page.labSolutions.title': 'Lab Solutions',
    'page.labSolutions.subtitle': "When a part is missing, you build it.",
    'page.park.title': 'Pseudomonas Park',
    'page.park.subtitle': 'The creative space.',
    'page.ideas.title': 'Ideas',
    'page.ideas.subtitle': 'Half-finished, wild, brewing.',
  },
} as const;

export type UIKey = keyof typeof ui[typeof defaultLang];
