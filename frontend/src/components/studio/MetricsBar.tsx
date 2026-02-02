import { useStudioStore } from '../../store/studioStore';

interface MetricProps {
    label: string;
    value: string | number | undefined;
    unit?: string;
    delta?: number;
    trend?: 'up' | 'down' | 'neutral';
}

function Metric({ label, value, unit, delta }: MetricProps) {
    const isPositiveDelta = delta && delta > 0;

    return (
        <div className="flex flex-col px-4 first:pl-0 border-r border-gray-100 last:border-0">
            <span className="text-[10px] font-black uppercase tracking-widest text-gray-400 mb-0.5">{label}</span>
            <div className="flex items-baseline gap-2">
                <span className="text-sm font-mono font-bold text-gray-800">
                    {value !== undefined ? value : '—'}{unit}
                </span>
                {delta !== undefined && delta !== 0 && (
                    <span className={`text-[10px] font-mono ${isPositiveDelta ? 'text-green-500' : 'text-red-500'}`}>
                        {isPositiveDelta ? '↑' : '↓'} {Math.abs(delta).toFixed(2)}
                    </span>
                )}
            </div>
        </div>
    );
}

export default function MetricsBar() {
    const { dashboard } = useStudioStore();

    if (!dashboard) return null;

    const { properties } = dashboard.baseline;
    const { property_delta, alerts } = dashboard;

    return (
        <div className="w-full bg-white border-b border-gray-200 px-6 py-3 flex items-center justify-between shadow-sm z-10">
            <div className="flex items-center">
                <Metric label="MW" value={properties.mw?.toFixed(1)} unit=" g/mol" />
                <Metric label="LogP" value={properties.logp?.toFixed(2)} delta={property_delta.logp} />
                <Metric label="TPSA" value={properties.tpsa?.toFixed(1)} unit=" Å²" delta={property_delta.tpsa} />
                <Metric label="QED" value={properties.qed?.toFixed(2)} delta={property_delta.qed} />
            </div>

            <div className="flex items-center gap-4">
                {alerts && alerts.length > 0 ? (
                    <div className="flex gap-2">
                        {alerts.slice(0, 3).map((alert, i) => (
                            <span
                                key={i}
                                className={`px-2 py-0.5 rounded text-[9px] font-black uppercase tracking-tighter ${alert.severity === 'high' ? 'bg-red-100 text-red-600' :
                                    alert.severity === 'medium' ? 'bg-amber-100 text-amber-600' :
                                        'bg-blue-100 text-blue-600'
                                    }`}
                            >
                                {alert.id}
                            </span>
                        ))}
                    </div>
                ) : (
                    <span className="text-[10px] font-black uppercase tracking-widest text-green-500 flex items-center gap-1">
                        <div className="w-1.5 h-1.5 rounded-full bg-green-500 animate-pulse" />
                        No Critical Alerts
                    </span>
                )}

                <div className="h-6 w-px bg-gray-100 mx-2" />

                <div className="flex flex-col items-end">
                    <span className="text-[8px] font-black uppercase tracking-[0.2em] text-gray-300">InChIKey</span>
                    <span className="text-[10px] font-mono text-gray-500">{dashboard.inchikey || 'PENDING...'}</span>
                </div>
            </div>
        </div>
    );
}
