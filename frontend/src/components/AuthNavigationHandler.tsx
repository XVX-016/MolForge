import { useEffect } from 'react';
import { useNavigate } from 'react-router-dom';
import { useAuthStore } from '../store/authStore';

export default function AuthNavigationHandler() {
  const navigate = useNavigate();
  const { user, intendedDestination } = useAuthStore();

  useEffect(() => {
    if (user && intendedDestination) {
      const destination = intendedDestination;
      useAuthStore.setState({ intendedDestination: null });
      navigate(destination, { replace: true });
    }
  }, [user, intendedDestination, navigate]);

  return null;
}








