/**
 * Animation utilities
 * 
 * Phase 15: Animation & UX Polish
 * 
 * Provides smooth transitions and easing functions.
 */

/**
 * Easing functions for animations
 */
export const easing = {
  linear: (t: number) => t,
  easeIn: (t: number) => t * t,
  easeOut: (t: number) => t * (2 - t),
  easeInOut: (t: number) => t < 0.5 ? 2 * t * t : -1 + (4 - 2 * t) * t,
  easeOutCubic: (t: number) => 1 - Math.pow(1 - t, 3),
  easeInCubic: (t: number) => t * t * t,
  easeInOutCubic: (t: number) => t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2,
  spring: (t: number) => {
    const c4 = (2 * Math.PI) / 3
    return t === 0 ? 0 : t === 1 ? 1 : Math.pow(2, -10 * t) * Math.sin((t * 10 - 0.75) * c4) + 1
  },
}

/**
 * Animate a value from start to end over duration
 */
export function animate(
  start: number,
  end: number,
  duration: number,
  easingFn: (t: number) => number = easing.easeOutCubic,
  onUpdate: (value: number) => void,
  onComplete?: () => void
): () => void {
  const startTime = performance.now()
  let animationFrame: number | null = null
  let cancelled = false

  const update = (currentTime: number) => {
    if (cancelled) return

    const elapsed = currentTime - startTime
    const progress = Math.min(elapsed / duration, 1)
    const eased = easingFn(progress)
    const value = start + (end - start) * eased

    onUpdate(value)

    if (progress < 1) {
      animationFrame = requestAnimationFrame(update)
    } else {
      if (onComplete) onComplete()
    }
  }

  animationFrame = requestAnimationFrame(update)

  // Return cancel function
  return () => {
    cancelled = true
    if (animationFrame !== null) {
      cancelAnimationFrame(animationFrame)
    }
  }
}

/**
 * Animate multiple values in parallel
 */
export function animateParallel(
  animations: Array<{
    start: number
    end: number
    duration: number
    easingFn?: (t: number) => number
    onUpdate: (value: number) => void
    onComplete?: () => void
  }>
): () => void {
  const cancels = animations.map(anim =>
    animate(
      anim.start,
      anim.end,
      anim.duration,
      anim.easingFn,
      anim.onUpdate,
      anim.onComplete
    )
  )

  return () => {
    cancels.forEach(cancel => cancel())
  }
}

/**
 * Debounce function with immediate option
 */
export function debounce<T extends (...args: any[]) => any>(
  func: T,
  wait: number,
  immediate: boolean = false
): (...args: Parameters<T>) => void {
  let timeout: NodeJS.Timeout | null = null

  return function executedFunction(...args: Parameters<T>) {
    const later = () => {
      timeout = null
      if (!immediate) func(...args)
    }

    const callNow = immediate && !timeout

    if (timeout) clearTimeout(timeout)
    timeout = setTimeout(later, wait)

    if (callNow) func(...args)
  }
}

/**
 * Throttle function
 */
export function throttle<T extends (...args: any[]) => any>(
  func: T,
  limit: number
): (...args: Parameters<T>) => void {
  let inThrottle: boolean

  return function executedFunction(...args: Parameters<T>) {
    if (!inThrottle) {
      func(...args)
      inThrottle = true
      setTimeout(() => (inThrottle = false), limit)
    }
  }
}

