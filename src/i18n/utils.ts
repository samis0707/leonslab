import { ui, defaultLang, type Lang, type UIKey } from './ui';

export function getLangFromUrl(url: URL): Lang {
  const [, maybeLang] = url.pathname.split('/');
  if (maybeLang in ui) return maybeLang as Lang;
  return defaultLang;
}

export function useTranslations(lang: Lang) {
  return function t(key: UIKey): string {
    return (ui[lang][key] ?? ui[defaultLang][key]) as string;
  };
}

/** Returns the URL path for the same page in the opposite language. */
export function getOppositeLangPath(url: URL): string {
  const lang = getLangFromUrl(url);
  const path = url.pathname;
  if (lang === 'de') {
    return '/en' + (path === '/' ? '/' : path);
  }
  const stripped = path.replace(/^\/en/, '') || '/';
  return stripped;
}

/** Returns a path prefixed with the current language (no prefix for default). */
export function localizedPath(lang: Lang, path: string): string {
  if (lang === defaultLang) return path;
  return `/${lang}${path}`;
}
